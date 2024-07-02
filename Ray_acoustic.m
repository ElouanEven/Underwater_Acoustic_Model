clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

SRC_z = -35;    % Source depth
REC_r = 4000;   % Receiver distance
REC_z = -40;    % Receiver depth
width = 20;     % Receiver width
f     = 1000;   % Frequency
theta_aperture = 40;  % Aperture angle (in degree)
num_rays = 50; % Number of rays to plot

surface_z = 0;  % Surface location
H = -200;       % Depth
W = REC_r*1.1;  % Maximum range

rho_air   = 1.225; % air density
rho_water = 1028;  % water density
rho_sand  = 1850;  % sand density

% Sound Speed Profile
depths   = [   0,  -10,  -20,  -35,  -65,  -80, -140, -200];
c_values = [1500, 1480, 1470, 1450, 1440, 1435, 1420, 1450]; %[1500, 1480, 1470, 1450, 1455, 1465, 1480, 1500]; %
c = @(z) interp1(depths, c_values, z, 'spline'); % Interpolated sound speed function
dt = 0.001;         % Time step
absorption = (3.3e-3 + (0.11*f^2)/(1+f^2) + (44*f^2)/(4100+f^2) +3e-4*f^2)*1e-3; % Absorption (dB per meter)

distance =   [ 0, W*3/4-1000, W*3/4-500, W*3/4, W];% [ 0, W/4, W/2, W*3/4, W];
H_values = H+[ 0,  50,   100,   50, 0];
H_d = @(r) interp1(distance, H_values, r, 'linear'); %'spline');


% Initialize intensity range for color mapping
log_intensity_min = inf;
log_intensity_max = -inf;

% Initialize rays vector
rays_r = cell(1, num_rays);
rays_z = cell(1, num_rays);
nbr_reflexion = zeros(1,num_rays);
rays_log_intensity = cell(1, num_rays);
theta_values = linspace(-theta_aperture/2, theta_aperture/2, num_rays);

% Plot of the SSP
Plot_SSP(c_values, c, depths, H, surface_z)

% Plot of the environment
Plot_Environment(H_values, H_d, distance, W, surface_z, SRC_z, REC_r, REC_z, width)

fprintf('----- Parameters : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Rays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Plot rays and compute intensities
for i = 1:num_rays          % For each rays
    theta = theta_values(i); % Initial angle at the source
    [rays_r{i}, rays_z{i}, rays_log_intensity{i}, nbr_reflexion(i)] = Ray_tracing(theta, SRC_z, H_d, c, W, absorption, dt);
    log_intensity_min = min(log_intensity_min, min(rays_log_intensity{i}));
    log_intensity_max = max(log_intensity_max, max(rays_log_intensity{i}));
end

% Rays arriving on the receiver
[ID_ray, delay_t, delay_i] = Cross_receiver(rays_r, rays_z, rays_log_intensity, REC_r, REC_z, width, num_rays, dt, nbr_reflexion);

% Rays arriving on the second receiver
[ID_ray_2, delay_t_2, delay_i_2] = Cross_receiver(rays_r, rays_z, rays_log_intensity, REC_r, REC_z-50, width, num_rays, dt, nbr_reflexion);

fprintf('----- Computing rays : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic


% % % For all rays without meaningful colors
% % Plot_All_Rays_Intensity(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, rays_log_intensity, log_intensity_min, log_intensity_max)


% For all rays without meaningful colors
Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z)


% For color intensity on the receiver
%Plot_Receiver_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, ID_ray, delay_i)

Plot_Receiver_Rays(H_d, W, surface_z, REC_r, REC_z-50, width, rays_r, rays_z, ID_ray_2, delay_i_2)

% Plot_Receiver_Delay(nbr_reflexion, ID_ray, delay_t, delay_i)

% Plot_Receiver_Delay(nbr_reflexion, ID_ray_2, delay_t_2, delay_i_2)


%% Cross validation with decreasing dt

ID_ray = ID_ray_2
nbr_CV = length(theta_values(ID_ray));

theta_values(ID_ray)
rays_r_improve = cell(1, nbr_CV);
rays_z_improve = cell(1, nbr_CV);
nbr_reflexion_improve = zeros(1,nbr_CV);
rays_log_intensity_improve = cell(1, nbr_CV);


for i = 1:nbr_CV
    theta = theta_values(ID_ray(i))
    [rays_r_improve{i}, rays_z_improve{i}, rays_log_intensity{i}, nbr_reflexion_improve(i)] = Ray_tracing(theta, SRC_z, H_d, c, W, absorption, 0.00001);
end
% Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r_improve, rays_z_improve)



figure;
hold on;
xlabel('Range (m)');
ylabel('Depth (m)');
title('Rac SSP all rays');
axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
% Plot surface and bottom
plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'r', 'LineWidth', 2);  % Bottom
plot([REC_r; REC_r], [REC_z-50; REC_z+width-50], 'r', 'LineWidth', 4);  % Receiver

% For all rays without color
for i = 1:nbr_CV
    r = rays_r{ID_ray(i)};
    z = rays_z{ID_ray(i)};
    line(r, z, 'Color', 'green'); % Plot all rays
end


for i = 1:nbr_CV
    r = rays_r_improve{i};
    z = rays_z_improve{i};
    line(r, z, 'Color', 'red'); % Plot all rays
end
hold off;

fprintf('----- Receiver : %.2f s -----\n', toc);
