clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Source and receiver position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SRC_z = -2;     % Source depth
REC_r = 10;     % Receiver distance %10
rec_z = -8;     % Receiver depth
width = 0.05; %0.1;    % Receiver width
space = 0.01; %0.05;    % Space between receivers
nbr_rec = 5;    % Number of receiver on the receiver array
REC_z = rec_z + (0:(nbr_rec-1)) * (width+space); % array of the z position of each receiver
theta_aperture = 140;  % Total aperture angle (in degree) %150
num_rays = 1500;  % Number of rays emmited from the source

surface_z = 0; % Surface location
H = -10;       % Depth
W = REC_r*1.1; % Maximum range  !!! Width depend to the placement of the receiver !!!

rho_air   = 1.225; % air density
rho_water = 1028;  % water density
rho_sand  = 1850;  % sand density

c_air = 340;    % Sound speed in air
c_water_mean = 1500; % Approximation of sound speed in water
c_sand = 1700;  % Sound speed in sand   !!! Page 39 Computational Ocean Acoustic !!!

% Sound Speed Profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depths   = [   0,   -1,   -2,   -3,   -4,   -5,   -6,   -7,   -8,   -9,  -10];
c_values = [1500, 1400, 1340, 1240, 1100, 1000, 1100, 1240, 1340, 1400, 1500]; % !!! Augmented values of real SSP !!!
c = @(z) interp1(depths, c_values, z, 'spline'); % Interpolated sound speed function

% Seabed altitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distance =   [ 0, W];
H_values = H+[ 0, 0];
H_d = @(r) interp1(distance, H_values, r, 'linear'); % Interpolated seabed altitude function

% Emitted signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amplitude_sig = 1; % Signal amplitude from the source
I_0 = (Amplitude_sig^2) / (2*rho_water*c_water_mean); % Source intensity
Fs = 2e5;     % Sampling frequency (Hz)
dt = 1/Fs;    % Time step
Fc = 1e4;     % Carrier frequency (Hz)
f = Fc*1e-3;  % f in kHz
absorption = (3.3e-3 + (0.11*f^2)/(1+f^2) + (44*f^2)/(4100+f^2) +3e-4*f^2)*1e-3; % Absorption (dB/m)    !!! Page 36 Computational Ocean Acoustic !!!

dur = 0.001;     % Signal duration (sec)

T = (0:(Fs*dur-1)) / Fs; % Time vector
sig = sin(2*pi*Fc*T); % Signal from the source
len_0 = Fc*2; % !!! Approximation of the total length vector !!!
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

% Plot of the SSP
Plot_SSP(c_values, c, depths, H, surface_z)

% Plot of the environment
Plot_Environment(H_values, H_d, distance, W, surface_z, SRC_z, REC_r, REC_z, width)

% Plot the emitted signal
figure;
hold on;
plot(T, sig, 'LineWidth', 2);
xlabel('Time (s)', 'LineWidth', 2);
ylabel('Amplitude', 'LineWidth', 2);
title(['Emitted signal with Fs = ', num2str(Fs), 'Hz and Fc = ', num2str(Fc), 'Hz'],  'LineWidth', 2);
hold off;


% Variable checkup

if REC_r < (2*Fc*((width+space)*nbr_rec)^2)/c(rec_z)
    disp("!!!!! Waves are considered as spherical !!!!!");
end

if Fs < 2*Fc
    disp("Temporal aliasing")
end

if width+space > c(rec_z)/(2*Fc)
    disp("Spatial aliasing")
end



fprintf('----- Parameters : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Rays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Plot rays and compute intensities
for i = 1:num_rays           % For each rays
    theta = theta_values(i); % Initial angle at the source
    [rays_r{i}, rays_z{i}, rays_intensity{i}, nbr_reflexion(i)] = Ray_tracing(theta, SRC_z, I_0, H_d, c, W, absorption, dt);
    intensity_min = min(intensity_min, min(rays_intensity{i}));
    intensity_max = max(intensity_max, max(rays_intensity{i}));
end

fprintf('----- Rays : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Cross_multiple_receivers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[ID_ray, delay_t, delay_i] = Cross_multiple_receivers(rays_r, rays_z, rays_intensity, REC_r, REC_z, width, dt, nbr_reflexion);

fprintf('----- Cross Receivers : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Time signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

SIG = zeros(len_0, nbr_rec);%cell(1, nbr_rec);
for k = 1:nbr_rec % For all receivers
    for i = 1:length(ID_ray{1,k}) % For all rays arriving on the receiver k
        i_timedelay = round(delay_t{1,k}(i)/dt);
        i_start = i_timedelay;
        i_end = i_timedelay+(length(sig));
        sig_delay = sig*sqrt(delay_i{1,k}(i) * 2*rho_water*c_water_mean); % Do we have to use the c_water_mean or the c(rec_z) ???
        SIG(i_start:i_end-1,k) = SIG(i_start:i_end-1,k) + sig_delay';
    end
end

fprintf('----- Time signal on receivers : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% For all rays without meaningful colors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

fprintf('----- Plot : %.2f s -----\n', toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Beamforming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Adding white noise to all signals

SNR = 50;

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot the original and filtered signals for comparison
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
title('Time signal with noise on receivers');
xlabel('Time (s)', 'LineWidth', 2);
ylabel('Amplitude', 'LineWidth', 2);
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Lfilt = 80;

N_phi = 180; % Phi resolution
PHI = 0:pi/N_phi:pi-pi/N_phi; % PHI is referenced to the vertical z-axis upward
% To obtain the angle of arrival on the receiver (theta) we should substract 90° 
c_iso = c(rec_z); % Sound speed considered iso around the receivers


%T_snap = 0:Lfilt:len_0-Lfilt;  % No overlap for now

Overlap = 0.6;          % Overlap
T_snap = 0:Lfilt*(1-Overlap):len_0-Lfilt; % Time for each start snapshot with overlap



% Frequency resolution
F = (0:Lfilt/2-1)*Fs/Lfilt;

ENERGY = zeros(length(T_snap),N_phi);
HannW = meshgrid(hann(Lfilt), zeros(nbr_rec,1)); % Hanning window to multiply to each x_snap sensor

for t_snap = 1:500 %length(T_snap)
    SIG_snap = SIG_noisy(T_snap(t_snap)+1:T_snap(t_snap)+Lfilt,:);
    SIG_snap_window = HannW' .* SIG_snap;       % Application of the hanning window on the snapshot signal
    FFT_x = fft(SIG_snap_window);             % Fourier transform
    Y = zeros(length(F), length(PHI));      % Initializing the output array
    
    for F_i = 1:length(F)
        for PHI_i = 1:length(PHI)
            w = zeros(nbr_rec,1);     % Initializing weights for all sensors 
            for rec_i = 1:nbr_rec
                u0 = [cos(PHI(PHI_i)); 0; 0];          % Direction normalised
                w(rec_i) = 1/nbr_rec * exp(-1i*2*pi/c_iso.*F(F_i)*[REC_z(rec_i)-REC_z(1) 0 0]*u0); % Weight equation
            end
            
            Y(F_i, PHI_i) = FFT_x(F_i,:) * w; % Output signal
        end
    end
    
    Y_abs = abs(Y);     % Absolute signal
    Energy = Y_abs.^2;  % Energy

    ENERGY(t_snap,:) = mean(Energy,1);
    
end



[T_mesh, PHIT_mesh] = meshgrid(T_snap/Fs, PHI*180/pi);

figure()
mesh(T_mesh, PHIT_mesh, ENERGY');
colormap(jet);
xlabel('Time (s)');
ylabel('Phi (deg)');
zlabel('Amplitude');
title('Energy (with BandPass filtering)');
set(gca, 'fontsize', 20, 'linewidth',2)
set(findall(gcf, 'Type', 'Line'),'LineWidth',2); 
xlim([0.006 0.03]);
colorbar;
