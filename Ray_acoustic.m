clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

Z_src = -50;    % Source depth
REC_r = 3500;   % Receiver distance
REC_z = -55;    % Receiver depth
width = 10;     % Receiver width
f     = 1000;   % Frequency
theta_aperture = 80;  % Aperture angle (in degree)
num_rays = 200; % Number of rays to plot

surface_z = 0;  % Surface location
H = -200;       % Depth
W = 4000;       % Maximum range

rho_air   = 1.225; % air density
rho_water = 1028;  % water density
rho_sand  = 1850;  % sand density

% Updated sound speed profile
depths   = [   0,  -10,  -20,  -35,  -65,  -80, -100, -200];
c_values = [1500, 1480, 1470, 1450, 1440, 1435, 1420, 1450];
c = @(z) interp1(depths, c_values, z, 'linear'); % Interpolated sound speed function
dt = 0.001;         % Time step
%absorption = 0.1;   % Absorption (dB per meter)
absorption = (3.3e-3 + (0.11*f^2)/(1+f^2) + (44*f^2)/(4100+f^2) +3e-4*f^2)*1e-3; 


% Initialize intensity range for color mapping
log_intensity_min = inf;
log_intensity_max = -inf;

% Initialize rays vector
rays_r = cell(1, num_rays);
rays_z = cell(1, num_rays);
nbr_reflexion = zeros(1,num_rays);
rays_log_intensity = cell(1, num_rays);
theta_values = linspace(-theta_aperture/2, theta_aperture/2, num_rays);

figure;
hold on;
plot(c_values, depths, 'o', 'LineWidth', 2);
plot(c(linspace(0, H, 200)), linspace(0, H, 200), 'LineWidth', 2); 
% Plot surface and bottom
plot([min(c(linspace(0, H, 200))), max(c(linspace(0, H, 200)))], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
plot([min(c(linspace(0, H, 200))), max(c(linspace(0, H, 200)))], [H, H], 'k', 'LineWidth', 2);  % Bottom
xlabel('Sound speed (m/s)', 'LineWidth', 2);
ylabel('Depth (m)', 'LineWidth', 2);
title('Sound Speed Profile', 'LineWidth', 2);
legend('Data', 'Interpolation', 'LineWidth', 2);
hold off;

fprintf('----- Parameters : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Ray %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Plot rays and compute intensities
for i = 1:num_rays          % For each rays
    theta = theta_values(i); % Initial angle at the source
    [rays_r{i}, rays_z{i}, rays_log_intensity{i}, nbr_reflexion(i)] = ray_trace(theta, Z_src, H, c, W, absorption, dt);
    log_intensity_min = min(log_intensity_min, min(rays_log_intensity{i}));
    log_intensity_max = max(log_intensity_max, max(rays_log_intensity{i}));
end

% Rays arriving on the receiver
[i_ray, delay_t, delay_i] = Cross_receiver(rays_r, rays_z, rays_log_intensity, REC_r , REC_z, width, num_rays, dt, nbr_reflexion);

fprintf('----- Computing rays : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Plot setup
figure;
hold on;
xlabel('Range (m)');
ylabel('Depth (m)');
title('Rac SSP receiver');
axis([0 W H surface_z]);

% Plot surface and bottom
plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
plot([0, W], [H, H], 'k', 'LineWidth', 2);  % Bottom
plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
% Color mapping based on log intensity
colormap('jet');
JET = jet(256);
for i = i_ray
    r = rays_r{i};
    z = rays_z{i};
    log_intensity = rays_log_intensity{i};
    color_idx = floor((log_intensity - log_intensity_min) / (log_intensity_max - log_intensity_min) * 255) + 1;
    color = JET(color_idx, :);  % Get color from colormap based on log intensity
    for j = 1:length(r)-1
        line([r(j), r(j+1)], [z(j), z(j+1)], 'Color', color(j,:));  % Apply color to line
    end
end


% Add color bar for intensity
caxis([log_intensity_min, log_intensity_max]);
colorbar;
ylabel(colorbar, 'Log Intensity');

hold off;


Delay_color = floor((nbr_reflexion(i_ray) - min(nbr_reflexion(i_ray))) / (max(nbr_reflexion(i_ray)) - min(nbr_reflexion(i_ray))) * 255) + 1;
figure;
hold on;
for i = 1:length(delay_t)
    plot(delay_t(i), delay_i(i), 'o', 'Color', JET(Delay_color(i), :), 'MarkerSize', 8, 'LineWidth', 2);
end
caxis([min(nbr_reflexion(i_ray)), max(nbr_reflexion(i_ray))]);
hold off;




fprintf('----- Receiver : %.2f s -----\n', toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, z, log_intensity, nbr_reflexion] = ray_trace(theta, Z_src, H, c, W, absorption, dt)
    
    % Initialize the ray path and intensity
    r = [    0,         c(Z_src)*cosd(theta)*dt];
    z = [Z_src, Z_src + c(Z_src)*sind(theta)*dt];
    intensity = [1, 10^(-absorption * c(z(end))*dt / 10)];
    count_ref = 0;
    
    while r(end) < W  % While the ray is in the rectangle
        
        % Actual conditions
        last_c    = c(z(end-1)); % Sound speed at last depth
        current_c = c(z(end  )); % Sound speed at current depth

        if abs((current_c/last_c)*cosd(theta)) < 1 % total transmission 
            theta = sign(theta)*acosd((current_c/last_c)*cosd(theta));
        else % total reflexion
            theta = -sign(theta)*acosd((last_c/current_c)*cosd(theta));
        end

        % Calculate increments in r and z
        dr = current_c * cosd(theta) * dt;
        dz = current_c * sind(theta) * dt;
        
        % Check for reflection
        if z(end) + dz > 0          % on surface
            r = [r, r(end) - z(end)/tand(theta),   r(end) + dr ];
            z = [z,                            0, -(z(end) + dz)];
            theta=-theta;
            intensity_ref = intensity(end) * 10^(-absorption * c(z(end))*dt / 10);
            intensity = [intensity, intensity_ref, intensity_ref];
            count_ref = count_ref + 1;

        elseif z(end) + dz < H      % on bottom
            r = [r, r(end) + (H-z(end))/tand(theta),       r(end) + dr ];
            z = [z,                                H, -(z(end)+dz) + 2*H];
            theta=-theta;
            intensity_ref = intensity(end) * 10^(-absorption * c(z(end))*dt / 10);
            intensity = [intensity, intensity_ref, intensity_ref];
            count_ref = count_ref + 1;

        else    % Nominal case
            r = [r, r(end) + dr];
            z = [z, z(end) + dz];
            intensity = [intensity, intensity(end) * 10^(-absorption * c(z(end))*dt / 10)];
        end
    end
    log_intensity = log10(intensity);
    nbr_reflexion = count_ref;
end


function [i_ray, delay_t, delay_i] = Cross_receiver(rays_r, rays_z, rays_log_intensity, REC_r , REC_z, width, num_rays, dt, nbr_reflexion)
    i_ray = [];
    delay_t = [];
    delay_i = [];
    for i = 1:num_rays
        for j = 1:length(rays_r{1,i})-1
            if (rays_r{1,i}(j) < REC_r) && (rays_r{1,i}(j+1) > REC_r)
                A = (rays_z{1,i}(j+1)-rays_z{1,i}(j))/(rays_r{1,i}(j+1)-rays_r{1,i}(j));
                B = rays_z{1,i}(j) - A*rays_r{1,i}(j);
                f_REC_r = A*REC_r + B;
                if (f_REC_r > REC_z) && (f_REC_r < REC_z+width)
                    i_ray = [i_ray, i];
                    delay_t = [delay_t, dt*(j - nbr_reflexion(i))];
                    delay_i = [delay_i, rays_log_intensity{1,i}(j)];
                end
            end
        end
    end
end
