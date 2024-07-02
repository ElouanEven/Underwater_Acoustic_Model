clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

Z_src = -35;    % Source depth
REC_r = 4000;   % Receiver distance
REC_z = -40;    % Receiver depth
width = 20;     % Receiver width
f     = 1000;   % Frequency
theta_aperture = 40;  % Aperture angle (in degree)
num_rays = 200; % Number of rays to plot

surface_z = 0;  % Surface location
H = -200;       % Depth
W = REC_r*1.1;  % Maximum range

rho_air   = 1.225; % air density
rho_water = 1028;  % water density
rho_sand  = 1850;  % sand density

% Sound Speed Profile
depths   = [   0,  -10,  -20,  -35,  -65,  -80, -140, -200];
c_values = [1500, 1480, 1470, 1450, 1455, 1465, 1480, 1500]; %[1500, 1480, 1470, 1450, 1440, 1435, 1420, 1450];
c = @(z) interp1(depths, c_values, z, 'spline'); % Interpolated sound speed function
dt = 0.003;         % Time step
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


figure;
hold on;
axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
plot(distance, H_values, 'o', 'LineWidth', 2);
plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'b', 'LineWidth', 2);
plot(0, Z_src, 'o', 'LineWidth', 5);                % Source
plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
% Plot surface and seabed
xlabel('Distance from source (m)', 'LineWidth', 2);
ylabel('Altitude (m)', 'LineWidth', 2);
title('Bathymetrie', 'LineWidth', 2);
legend('Data', 'Interpolation', 'Source', 'Receiver', 'LineWidth', 2);
hold off;




fprintf('----- Parameters : %.2f s -----\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Rays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Plot rays and compute intensities
for i = 1:num_rays          % For each rays
    theta = theta_values(i); % Initial angle at the source
    [rays_r{i}, rays_z{i}, rays_log_intensity{i}, nbr_reflexion(i)] = ray_trace(theta, Z_src, H_d, c, W, absorption, dt);
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

% For all rays without color

% Plot setup
figure;
hold on;
xlabel('Range (m)');
ylabel('Depth (m)');
title('Rac SSP all rays');
axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
% Plot surface and bottom
plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'r', 'LineWidth', 2);  % Bottom
plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
% Color mapping based on log intensity
colormap('gray');
JET = gray(256);

% % % For all rays with color intensity
% % for i = 1:num_rays
% %     r = rays_r{i};
% %     z = rays_z{i};
% %     log_intensity = rays_log_intensity{i};
% % 
% %    for j = 1:length(r)-1
% %         color_idx = floor((log_intensity(j) - log_intensity_min) / (log_intensity_max - log_intensity_min) * 255) + 1;
% %         color = jet(256);
% %         color = color(color_idx, :);  % Get color from colormap based on log intensity
% %         line([r(j), r(j+1)], [z(j), z(j+1)], 'Color', color);  % Apply color to line
% %     end
% % end

% For all rays without color
for i = 1:num_rays
    r = rays_r{i};
    z = rays_z{i};
    log_intensity = rays_log_intensity{i};
    line(r, z);%, 'Color', color);  % Apply color to line
end

% Add color bar for intensity
caxis([min(delay_i), max(delay_i)]);
colorbar;
ylabel(colorbar, 'Log Intensity');

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For color intensity on the receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot setup
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
% Color mapping based on log intensity
colormap('gray');
JET = gray(256);

j=0;
for i = i_ray
    r = rays_r{i};
    z = rays_z{i};
    log_intensity = rays_log_intensity{i};
    j=j+1;
    color_idx = 256 - floor((delay_i(j) - min(delay_i)) / (max(delay_i) - min(delay_i)) * 255);
    color = JET(color_idx, :);  % Get color from colormap based on log intensity
    line(r, z, 'Color', color);  % Apply color to line
end



% Add color bar for intensity
caxis([min(delay_i), max(delay_i)]);
colorbar;
ylabel(colorbar, 'Log Intensity');

hold off;




if max(nbr_reflexion(i_ray)) == min(nbr_reflexion(i_ray))
    Delay_color = ones(1,length(delay_t));
else
    Delay_color = floor((nbr_reflexion(i_ray) - min(nbr_reflexion(i_ray))) / (max(nbr_reflexion(i_ray)) - min(nbr_reflexion(i_ray))) * 255) + 1;
end

figure;
hold on;
for i = 1:length(delay_t)
    plot(delay_t(i), delay_i(i), 'o', 'Color', JET(Delay_color(i), :), 'MarkerSize', 8, 'LineWidth', 2);
end
hold off;




fprintf('----- Receiver : %.2f s -----\n', toc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, z, log_intensity, nbr_reflexion] = ray_trace(theta, Z_src, H_d, c, W, absorption, dt)
    MAX_STEP = floor(W/(c(Z_src)*dt) * 1.5); % Estimation of the max number of increments with theta_aperture < 90
    r = zeros(1,MAX_STEP);
    z = zeros(1,MAX_STEP);
    intensity = zeros(1,MAX_STEP);

    % Initialize the ray path and intensity
    r(1:2) = [    0,         c(Z_src)*cosd(theta)*dt];
    z(1:2) = [Z_src, Z_src + c(Z_src)*sind(theta)*dt];
    intensity(1:2) = [1, 10^(-absorption * c(z(end))*dt / 10)];
    count_ref = 0;
    i = 2;
    while (r(end) < W && i < MAX_STEP)  % While the ray is in the rectangle
        
        % Actual conditions
        last_c    = c(z(i-1)); % Sound speed at last depth
        current_c = c(z(i  )); % Sound speed at current depth

        if abs((current_c/last_c)*cosd(theta)) < 1 % total transmission 
            theta = sign(theta)*acosd((current_c/last_c)*cosd(theta));
        else % total reflexion
            theta = -sign(theta)*acosd((last_c/current_c)*cosd(theta));
        end
        
        % Calculate increments in r and z
        dr = current_c * cosd(theta) * dt;
        dz = current_c * sind(theta) * dt;
        
        % Check for reflection
        if z(i) + dz > 0          % on surface
            r(i+1:i+2) = [r(i) - z(i)/tand(theta), r(i) + dr ];
            z(i+1:i+2) = [0, -(z(i) + dz)];
            theta=-theta;
            intensity_ref = intensity(i) * 10^(-absorption * c(z(i))*dt / 10);
            intensity(i+1:i+2) = [intensity_ref, intensity_ref];
            count_ref = count_ref + 1;
            i=i+2;

        elseif z(i) + dz < H_d(r(i) + dr)      % on bottom
            
            theta_bot = atand(H_d(r(i)+1)-H_d(r(i))); % Bottom slop (calcul with 1m)
            r_bot = r(i) + (z(i)-H_d(r(i)))/(tand(abs(theta))+tand(theta_bot));
            z_bot = H_d(r_bot);
            dist_ground = sqrt((r(i)+dr - r_bot)^2 + (z(i)+dz - z_bot)^2);
            theta = abs(theta)+2*theta_bot;

            if theta > 90 % No rays rackward
                r = r(1:i);
                z = z(1:i);
                intensity = intensity(1:i);
                log_intensity = log10(intensity);
                nbr_reflexion = count_ref;
                return
            end

            r_ref = r_bot + dist_ground*cosd(theta);
            z_ref = z_bot + dist_ground*sind(theta);
            r(i+1:i+2) = [r_bot, r_ref];
            z(i+1:i+2) = [z_bot, z_ref];
            
            intensity_ref = intensity(i) * 10^(-absorption * c(z(i))*dt / 10);
            intensity(i+1:i+2) = [intensity_ref, intensity_ref];
            count_ref = count_ref + 1;
            i=i+2;

        else    % Nominal case
            r(i+1) = r(i) + dr;
            z(i+1) = z(i) + dz;
            intensity(i+1) = intensity(i) * 10^(-absorption * c(z(i))*dt / 10);
            i=i+1;
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
