clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Source and receiver position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SRC_z = -2;     % Source depth
REC_r = 50;     % Receiver distance
rec_z = -8;     % Receiver depth
width = 0.3;    % Receiver width
space = 0.1;    % Space between receivers
nbr_rec = 3;    % Number of receiver on the receiver array
REC_z = rec_z + (0:(nbr_rec-1)) * (width+space); % array of the z position of each receiver
theta_aperture = 150;  % Total aperture angle (in degree)
num_rays = 200;  % Number of rays emmited from the source

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
Fs = 1e5;     % Sampling frequency (Hz)
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

% For all rays without meaningful colors
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z)

fprintf('----- Rays : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Cross_multiple_receivers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

[ID_ray, delay_t, delay_i] = Cross_multiple_receivers(rays_r, rays_z, rays_intensity, REC_r, REC_z, width, dt, nbr_reflexion);

fprintf('----- Receiver : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Plot_Rays_on_all_receivers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
xlabel('Range (m)');
ylabel('Depth (m)');
title('Rac SSP receiver');
axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
% Plot surface and bottom
plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'r', 'LineWidth', 2);  % Bottom
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
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Delays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JET = jet(256);

figure;
hold on;
for k = 1:nbr_rec
    Delay_color = floor((k / nbr_rec) * 255) + 1;
    plot(delay_t{1,k}, delay_i{1,k}, 'o', 'Color', JET(Delay_color, :), 'MarkerSize', 8, 'LineWidth', 2);
end

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Time signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SIG = cell(1, nbr_rec);

for k = 1:nbr_rec
    SIG{1,k} = zeros(1,len_0);

    for i = 1:length(ID_ray{1,k})
        i_timedelay = round(delay_t{1,k}(i)/dt);
        i_start = i_timedelay;
        i_end = i_timedelay+(length(sig));
        sig_delay = sig*sqrt(delay_i{1,1}(i) * 2*rho_water*c_water_mean);
        SIG{1,k}(i_start:i_end-1) = SIG{1,k}(i_start:i_end-1) + sig_delay;
    end
end

figure;
hold on
for k = 1:nbr_rec
    Delay_color = floor((k / nbr_rec) * 255) + 1;
    plot(T_0, SIG{1,k}, 'Color', JET(Delay_color, :));
end
hold off

fprintf('----- Computing rays : %.2f s -----\n', toc);
