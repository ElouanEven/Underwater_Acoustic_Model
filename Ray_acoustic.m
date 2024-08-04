clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Source and receiver position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SRC_z = -2;     % Source depth
REC_r = 10;     % Receiver distance
rec_z = -8;     % Receiver depth
width = 0.04;   % Receiver width
space = 0.005;  % Space between receivers
nbr_rec = 10;   % Number of receiver on the receiver array
REC_z = rec_z + (0:(nbr_rec-1)) * (width+space); % array of the z position of each receiver
REC_z_centered = mean(REC_z); % Center of the receiver
theta_aperture = 100;  % Total aperture angle (in degree)
num_rays = 1000;  % Number of rays emmited from the source

% Environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surface_z = 0; % Surface location
H = -10;       % Depth
W = REC_r*1.1; % Maximum range

rho_air   = 1.225; % air density
rho_water = 1028;  % water density
rho_sand  = 1850;  % sand density   !!! Defining real rho !!!

c_air = 340;    % Sound speed in air
c_water_mean = 1500; % Approximation of sound speed in water
c_sand = 1700;  % Sound speed in sand   !!! Page 39 Computational Ocean Acoustic !!!

% Sound Speed Profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depths   = [   0,   -1,   -2,   -3,   -4,   -5,   -6,   -7,   -8,   -9,  -10];
c_values =  [1500, 1400, 1340, 1240, 1100, 1000, 1100, 1240, 1300, 1350, 1400]; %[1500, 1480, 1470, 1450, 1440, 1435, 1420, 1410, 1405, 1400, 1400]; %[1500, 1400, 1340, 1240, 1100, 1000, 1100, 1240, 1340, 1400, 1500]; %[1500, 1470, 1400, 1350, 1280, 1200, 1100, 1000, 970, 960, 950]; % !!! Augmented values of real SSP !!!
c = @(z) interp1(depths, c_values, z, 'spline'); % Interpolated sound speed function

% Approximation of seabed reflexion in dB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_no_reflexion = atand(sqrt((1-(c_sand/c_water_mean)^2)/(((rho_sand*c_sand)/(rho_water*c_water_mean)^2)-1))); % Maximum value of theta before all energy is transmitted to the bottom !!! c_water_mean or c(H_d(r(i)) ??? !!!
Angle = [ 0, -acosd(c_water_mean/c_sand), -90];
R_values = [0, 0, 10]; 
R_coef = @(z) interp1(Angle, R_values, z, 'linear');

% Seabed altitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distance =   [ 0, 2, 3, 5, 7, W];
H_values = H+[ 0, 0, 0, 0, 0, 0];
H_d = @(r) interp1(distance, H_values, r, 'spline'); %'linear'); % Interpolated seabed altitude function

% Emitted signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amplitude_sig = 1; % Signal amplitude from the source
I_0 = (Amplitude_sig^2) / (2*rho_water*c_water_mean); % Source intensity   !!! Which value of c do we have to use ??? !!!
Fs = 5e5;     % Sampling frequency (Hz)

dt = 1/1e4;   % Time step     !!! Not linked to the Sampling Frequency anymore !!!
Fc = 1e4;     % Carrier frequency (Hz)
f = Fc*1e-3;  % f in kHz
absorption = (3.3e-3 + (0.11*f^2)/(1+f^2) + (44*f^2)/(4100+f^2) +3e-4*f^2)*1e-3; % Absorption (dB/m)    !!! Page 36 Computational Ocean Acoustic !!!
              % Compute absorption : http://resource.npl.co.uk/acoustics/techguides/seaabsorption/

dur = 0.0005;  % Signal duration (sec)

T = (0:(Fs*dur-1)) / Fs; % Time vector
sig = sin(2*pi*Fc*T); % Signal from the source
len_0 = floor(2*W/min(c(-10:0))*Fs); % !!! Approximation of the total time length vector !!!
T_0 = (0:(len_0)-1) / Fs;

% Initialize intensity range for color mapping
intensity_min = inf;
intensity_max = 0;

% Initialize rays vector
rays_r = cell(1, num_rays);
rays_z = cell(1, num_rays);
nbr_reflexion = zeros(1,num_rays);
rays_intensity = cell(1, num_rays);
theta_values = linspace(-theta_aperture/2, theta_aperture/2, num_rays);

% Check for valid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if REC_r < (2*Fc*((width+space)*nbr_rec)^2)/c(rec_z)
    disp("!!!!! Waves are considered spherical !!!!!");
end

if Fs < 2*Fc
    disp("Temporal aliasing")
end

if width+space > c(rec_z)/(2*Fc)
    disp("Spatial aliasing")
end

if dur < 2/Fc
    disp("Signal duration too short")
end

fprintf('----- Parameters : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Rays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Plot rays and compute intensities
for i = 1:num_rays           % For each rays
    theta = theta_values(i); % Initial angle at the source
    [rays_r{i}, rays_z{i}, rays_intensity{i}, nbr_reflexion(i)] = Ray_tracing(theta, SRC_z, I_0, H_d, c, W, absorption, dt, c_sand, rho_sand, rho_water, R_coef);
    intensity_min = min(intensity_min, min(rays_intensity{i}));
    intensity_max = max(intensity_max, max(rays_intensity{i}));
end

fprintf('----- Rays : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Cross_multiple_receivers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[ID_ray, delay_t, delay_i] = Cross_multiple_receivers(rays_r, rays_z, rays_intensity, REC_r, REC_z, width, dt, nbr_reflexion, c(rec_z), absorption);

fprintf('----- Cross Receivers : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Time signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

SIG = zeros(len_0, nbr_rec);
for k = 1:nbr_rec % For all receivers
    for i = 1:length(ID_ray{1,k}) % For all rays arriving on the receiver k
        i_timedelay = round(delay_t{1,k}(i)*Fs);
        i_start = i_timedelay;
        i_end = i_timedelay+(length(sig));
        sig_delay = sig*sqrt(delay_i{1,k}(i) * 2*rho_water*c_water_mean); % Do we have to use the c_water_mean or the c(rec_z) ???
        if i_end < len_0
            SIG(i_start:i_end-1,k) = SIG(i_start:i_end-1,k) + sig_delay';
        else
            fprintf("Issue with ray : %.5f \n", delay_t{1,k}(i)); % The ray didn't reach the receiver before the last step ==> Increase the value of MAX_STEP
        end
    end
end

% Adding white noise to all signals
SNR = 0;
Amplitude_noise = Amplitude_sig /(10^(SNR/20)); % noise power at each sensor
NOISE = sqrt(Amplitude_noise)*randn(len_0,nbr_rec);
SIG_noisy = SIG + NOISE;

% Define the passband frequencies
f_low = 8000;  % Lower cutoff frequency in Hz
f_high = 12000; % Upper cutoff frequency in Hz

% Normalize the frequencies by the Nyquist frequency (Fs/2)
nyquist = Fs / 2;
Wn = [f_low f_high] / nyquist;

% Design the 4th-order Butterworth filter
[b, a] = butter(4, Wn, 'bandpass');

% Apply the filter to the signal using filtfilt for zero-phase filtering
filtered_SIG = filtfilt(b, a, SIG_noisy);

fprintf('----- Time signal on receivers with noise : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Beamforming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

Lfilt = 100; % Frequency resolution
N_phi = 180; % Phi resolution
PHI = 0:180/N_phi:180-180/N_phi; % PHI is referenced to the vertical z-axis upward
% To obtain the angle of arrival on the receiver (theta) we should substract 90°
c_iso = c(rec_z); % Sound speed considered iso around the receivers

Overlap = 0.3;          % Overlap
T_snap = 0:Lfilt*(1-Overlap):len_0-Lfilt; % Time for each start snapshot with overlap
[T_mesh, PHIT_mesh] = meshgrid(T_snap/Fs, PHI);
threshold = 10; % Min amplitude to consider a direction

% Frequency resolution
F = (0:Lfilt/2-1)*Fs/Lfilt;

[ENERGY, list_angle_detected, list_time_detected] = Beamformer_DaS(filtered_SIG, REC_z, Lfilt, F, PHI, c_iso, T_0, T_snap, threshold);

fprintf('----- Beamformer : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finding Source (ray tracing backward) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

rays_r_det = cell(1, length(list_angle_detected));
rays_z_det = cell(1, length(list_angle_detected));

for j=1:length(list_angle_detected)
    theta_detec = list_angle_detected(j);
    time_detec  = list_time_detected(j);
    
    [r, z] = Finding_Source(theta_detec, time_detec, REC_r, REC_z_centered, H_d, c, W, dt);
    
    rays_r_det{j} = r;
    rays_z_det{j} = z;
end

fprintf('----- Finding Source : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Plot of the SSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_SSP(c_values, c, depths, H, surface_z)

% Plot of the environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_Environment(H_values, H_d, distance, W, surface_z, SRC_z, REC_r, REC_z, width)

% Plot the emitted signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
plot(T, sig, 'LineWidth', 2);
xlabel('Time (s)', 'LineWidth', 2);
ylabel('Amplitude', 'LineWidth', 2);
title(['Emitted signal with Fs = ', num2str(Fs), 'Hz and Fc = ', num2str(Fc), 'Hz'],  'LineWidth', 2);
hold off;

% % % For all rays with intensity gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot_All_Rays_Intensity(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, rays_intensity, intensity_min, intensity_max)

% For all rays with meaningless colors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z)

% For all rays touching the receivers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_Rays_on_all_receivers(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, ID_ray)

% Plot rescale on receivers
figure;
hold on;
xlabel('Range (m)');
ylabel('Depth (m)');
title('Rac SSP receiver (zoomed)');
plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
% Color mapping based on intensity
colormap('jet');
JET = jet(256);
for k = 1:nbr_rec
    color = floor((k / nbr_rec) * 255) + 1;
    for i = ID_ray{1,k}
        line(rays_r{1,i}, rays_z{1,i}, 'Color', JET(color, :));
    end
end
H_imag = REC_z(end)+width+0.1 - REC_z(1)-0.1;
xlim([REC_r-H_imag REC_r+H_imag]);
ylim([REC_z(1)-0.1 REC_z(end)+width+0.1]);
hold off;

% Delays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JET = jet(256);
figure;
hold on;
for k = 1:nbr_rec
    Delay_color = floor((k / nbr_rec) * 255) + 1;
    plot(delay_t{1,k}, delay_i{1,k}, 'o', 'Color', JET(Delay_color, :), 'MarkerSize', 8, 'LineWidth', 2);
end
title('Intensity function of time delay on receivers');
xlabel('Time delay (s)', 'LineWidth', 2);
ylabel('Intensity (W/m^2)', 'LineWidth', 2);
hold off;

% Time signal on receivers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on
for k = 1:nbr_rec
    Delay_color = floor((k / nbr_rec) * 255) + 1;
    plot(T_0, SIG(:,k), 'Color', JET(Delay_color, :));
end
title('Time signal on receivers');
xlabel('Time (s)', 'LineWidth', 2);
ylabel('Amplitude', 'LineWidth', 2);
hold off

% Plot the original and filtered signals for comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2, 1, 1);
hold on
for k = 1:nbr_rec
    Delay_color = floor((k / nbr_rec) * 255) + 1;
    plot(T_0, SIG_noisy(:,k), 'Color', JET(Delay_color, :));
end
title('Time signal with noise on receivers');
xlabel('Time (s)', 'LineWidth', 2);
ylabel('Amplitude', 'LineWidth', 2);
hold off

subplot(2, 1, 2);
hold on
for k = 1:nbr_rec
    Delay_color = floor((k / nbr_rec) * 255) + 1;
    plot(T_0, filtered_SIG(:,k), 'Color', JET(Delay_color, :));
end
title('Time signal with noise filtered');
xlabel('Time (s)', 'LineWidth', 2);
ylabel('Amplitude', 'LineWidth', 2);
hold off

% DaS Beamformer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
mesh(T_mesh, PHIT_mesh, ENERGY');
colormap(jet);
xlabel('Time (s)');
ylabel('Phi (deg)');
zlabel('Amplitude');
title('Energy (with BandPass filtering)');
set(gca, 'fontsize', 20, 'linewidth',2)
set(findall(gcf, 'Type', 'Line'),'LineWidth',2);
colorbar;

% For all rays from the receiver with meaningless colors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r_det, rays_z_det)

% Link between both direction (from the source and the receiver)
Plot_Difference_estimation_reality(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, ID_ray, rays_r_det, rays_z_det)

fprintf('----- Plot : %.2f s -----\n', toc);
