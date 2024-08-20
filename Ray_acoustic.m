clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Source and receiver position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SRC_z = -5;     % Source depth
REC_r = 200;    % Receiver distance
rec_z = -40;    % Receiver depth
width = 0.04;   % Receiver width
space = 0.005;  % Space between receivers
nbr_rec = 20;   % Number of receiver on the receiver array
REC_z = rec_z + (0:(nbr_rec-1)) * (width+space); % array of the z position of each receiver
REC_z_centered = mean(REC_z); % Center of the receiver
theta_aperture = 90;  % Total aperture angle (in degree)
num_rays = 180;  % Number of rays emmited from the source

ZONE_width = (REC_z(end)-REC_z(1))+ 30*width; % Width of the first estimation of angle research (blue zone)
theta_resolution = theta_aperture/num_rays; % Resolution angle of the emmited rays from the source
coef_improvment_Dthteta = 40; % Coefficient of angular resolution improvment for ray emmission

% Environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surface_z = 0; % Surface location
W = REC_r*1.1; % Maximum range

rho_air   = 1.225; % air density
rho_water = 1028;  % water density
rho_sand  = 1850;  % sand density   !!! Defining real rho !!!

c_air = 340;    % Sound speed in air
c_water_mean = 1500; % Approximation of sound speed in water
c_sand = 1700;  % Sound speed in sand   !!! Page 39 Computational Ocean Acoustic !!!

% Real Sound Speed Profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = -50;
depths   = [   0,   -5,  -10,  -15,  -20,  -25,  -30,  -35,  -40,  -45,  -50];
c_values = [1500, 1480, 1470, 1450, 1440, 1435, 1420, 1410, 1405, 1400, 1400]; % [1500, 1400, 1340, 1240, 1100, 1000, 1100, 1240, 1340, 1400, 1500]; %

c = @(z) interp1(depths, c_values, z, 'spline'); % Interpolated sound speed function

% Approximation of seabed reflexion in dB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_no_reflexion = atand(sqrt((1-(c_sand/c_water_mean)^2)/(((rho_sand*c_sand)/(rho_water*c_water_mean)^2)-1))); % Maximum value of theta before all energy is transmitted to the bottom !!! c_water_mean or c(H_d(r(i)) ??? !!!
Angle = [ 0, -acosd(c_water_mean/c_sand), -90];
R_values = [0, 0, 10];
R_coef = @(z) interp1(Angle, R_values, z, 'linear');

% Seabed altitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distance = linspace(0, W, 6);
H_values = H+[ 0, 0, 0, 0, 0, 0];
H_d = @(r) interp1(distance, H_values, r, 'spline'); %'linear'); % Interpolated seabed altitude function

% Emitted signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amplitude_sig = 1; % Signal amplitude from the source
I_0 = (Amplitude_sig^2) / (2*rho_water*c_water_mean); % Source intensity   !!! Which value of c do we have to use ??? !!!
Fs = 5e5;    % Sampling frequency (Hz)

dt = 1/1e4;  % Time step     !!! Not linked to the Sampling Frequency anymore !!!
Fc = 1e4;    % Carrier frequency (Hz)
f = Fc*1e-3; % f in kHz
absorption = (3.3e-3 + (0.11*f^2)/(1+f^2) + (44*f^2)/(4100+f^2) +3e-4*f^2)*1e-3; % Absorption (dB/m)    !!! Page 36 Computational Ocean Acoustic !!!
             % Compute absorption : http://resource.npl.co.uk/acoustics/techguides/seaabsorption/

dur = 0.0005; % Signal duration (sec)
T = (0:(Fs*dur-1)) / Fs; % Time vector
sig = sin(2*pi*Fc*T); % Signal from the source
len_0 = floor(4*W/min(c(-10:0))*Fs); % !!! Approximation of the total time length vector !!!
T_0 = (0:(len_0)-1) / Fs;

% Define the passband frequencies
f_low = 8000;   % Lower cutoff frequency in Hz
f_high = 12000; % Upper cutoff frequency in Hz
nyquist = Fs / 2; % Nyquist frequency

% Initialize rays vector
rays_r = cell(1, num_rays+1);
rays_z = cell(1, num_rays+1);
nbr_reflexion = zeros(1,num_rays+1);
rays_intensity = cell(1, num_rays+1);
theta_values = linspace(-theta_aperture/2, theta_aperture/2, num_rays+1); %List of ray angle from the source

% Beamformer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lfilt = 200; % Frequency resolution
F = 0:Fs/Lfilt:Fs/2-Fs/(2*Lfilt); % Frequency vector
N_phi = 360; % Phi resolution
PHI = 0:180/N_phi:180-180/N_phi; % PHI is referenced to the vertical z-axis upward
% To obtain the angle of arrival on the receiver (theta) we should substract 90°
c_iso = c(rec_z); % Sound speed considered iso around the receivers

Overlap = 0.7; % Overlap
T_snap = 0:round(Lfilt*(1-Overlap)):len_0-Lfilt; % Time for each start snapshot with overlap
threshold = 2; % Min amplitude to consider a direction
%
%
%


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Rays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Plot rays and compute intensities
for i = 1:num_rays           % For each rays
    theta = theta_values(i); % Initial angle at the source
    [rays_r{i}, rays_z{i}, rays_intensity{i}, nbr_reflexion(i)] = Ray_tracing(theta, SRC_z, I_0, H_d, c, W, absorption, dt, R_coef);
end

fprintf('----- Rays : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Cross_zone %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ID_ray_zone] = Cross_zone(rays_r, rays_z, REC_r, REC_z_centered, ZONE_width); % ID of the rays crossing the zone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Rays_arround_zone %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[theta_values_zone] = Theta_Improved(ID_ray_zone, theta_values, theta_resolution, coef_improvment_Dthteta);

rays_r_zone = cell(1, length(theta_values_zone));
rays_z_zone = cell(1, length(theta_values_zone));
rays_intensity_zone = cell(1, length(theta_values_zone));
nbr_reflexion_zone = zeros(1, length(theta_values_zone));

% Plot rays and compute intensities
for i = 1:length(theta_values_zone) % For each rays
    theta = theta_values_zone(i); % Initial angle at the source
    [rays_r_zone{i}, rays_z_zone{i}, rays_intensity_zone{i}, nbr_reflexion_zone(i)] = Ray_tracing(theta, SRC_z, I_0, H_d, c, W, absorption, dt, R_coef);
end

fprintf('----- Rays arround zone : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Cross_multiple_receivers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ID_ray, delay_t, delay_i] = Cross_multiple_receivers(rays_r_zone, rays_z_zone, rays_intensity_zone, REC_r, REC_z, width, dt, nbr_reflexion_zone, c(rec_z), absorption);
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r_zone, rays_z_zone, ZONE_width)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Time signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[SIG] = Create_SIG(nbr_rec, len_0, ID_ray, sig, Fs, delay_t, delay_i, rho_water, c_water_mean);

% Adding white noise to all signals
SNR = 0;
Amplitude_noise = Amplitude_sig /(10^(SNR/20)); % noise power at each sensor
NOISE = sqrt(Amplitude_noise)*randn(len_0,nbr_rec);
SIG_noisy = SIG + NOISE;

Wn = [f_low f_high] / nyquist; % Normalize the frequencies
[b, a] = butter(4, Wn, 'bandpass'); % Design the 4th-order Butterworth filter
filtered_SIG = filtfilt(b, a, SIG_noisy); % Apply the filter to the signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Beamforming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[ENERGY, list_angle_detected, list_time_detected] = Beamformer_DaS(filtered_SIG, REC_z, Lfilt, F, PHI, c_iso, T_0, T_snap);

fprintf('----- Beamformer : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finding Source (ray tracing backward) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

rays_r_det = cell(1, length(list_angle_detected));
rays_z_det = cell(1, length(list_angle_detected));

for j=1:length(list_angle_detected)
    theta_detec = list_angle_detected(j);
    
    [r, z] = Finding_Source(theta_detec, REC_r, REC_z_centered, H_d, c, W, dt);
    
    rays_r_det{j} = r;
    rays_z_det{j} = z;
end

[r_inter, z_inter, t_match, Best_time, R_SRC, Z_SRC] = Find_Intersection(rays_r_det, rays_z_det, dt, list_time_detected, list_angle_detected);

fprintf('----- Finding Source : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% % % Plot of the SSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot_SSP(c_values, c, depths, H, surface_z)
% % 
% % % Plot of the environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot_Environment(H_values, H_d, distance, W, surface_z, SRC_z, REC_r, REC_z, width)
% % 
% % % Plot the emitted signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure;
% % hold on;
% % plot(T, sig, 'LineWidth', 2);
% % xlabel('Time (s)', 'LineWidth', 2);
% % ylabel('Amplitude', 'LineWidth', 2);
% % title(['Emitted signal with Fs = ', num2str(Fs), 'Hz and Fc = ', num2str(Fc), 'Hz'],  'LineWidth', 2);
% % hold off;
% % 
% % % For all rays with intensity gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot_All_Rays_Intensity(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, rays_intensity, intensity_min, intensity_max)

% For all rays with meaningless colors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, ZONE_width)
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r(ID_ray_zone), rays_z(ID_ray_zone), ZONE_width)
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r_zone, rays_z_zone, ZONE_width)

% For all rays touching the receivers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_Rays_on_all_receivers(H_d, W, surface_z, REC_r, REC_z, width, rays_r_zone, rays_z_zone, ID_ray)
Plot_Rays_zoomed_receivers(REC_r, REC_z, width, rays_r_zone, rays_z_zone, ID_ray)

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
[T_mesh, PHIT_mesh] = meshgrid(T_snap(1:length(ENERGY(:,1)))/Fs, PHI);

figure()
mesh(T_mesh, PHIT_mesh, ENERGY');
colormap(jet);
xlabel('Time (s)');
ylabel('Phi (deg)');
zlabel('Amplitude');
title('Energy (with BandPass filtering)');
set(gca, 'fontsize', 20, 'linewidth',2)
set(findall(gcf, 'Type', 'Line'),'LineWidth',2);
set(gca, 'ZScale', 'log');
colorbar;

% For all rays from the receiver with meaningless colors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r_det, rays_z_det, ZONE_width)

% Link between both direction (from the source and the receiver)
Plot_Difference_estimation_reality(H_d, W, surface_z, REC_r, REC_z, width, rays_r_zone, rays_z_zone, ID_ray, rays_r_det, rays_z_det, r_inter, z_inter)

ERROR = sqrt((0-R_SRC)^2+(SRC_z-Z_SRC));
fprintf('The final estimation error is %.3f m at the point {%.2f,%.2f} \n', ERROR, r_inter(1), z_inter(1));
fprintf('# --- {  r ,  z  } ;   precision --- \n');
fprintf('1  :  {%.2f,%.2f} ;   %.2e s \n', r_inter(1), z_inter(1), t_match(1));
fprintf('2  :  {%.2f,%.2f} ;   %.2e s \n', r_inter(2), z_inter(2), t_match(2));
fprintf('3  :  {%.2f,%.2f} ;   %.2e s \n', r_inter(3), z_inter(3), t_match(3));
fprintf('4  :  {%.2f,%.2f} ;   %.2e s \n', r_inter(4), z_inter(4), t_match(4));
