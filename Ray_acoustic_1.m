clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%atan(sqrt((1-(current_c/last_c)^2)/(((rho_air*current_c)/(rho_water*last_c))^2-1)))
Z_src = -50;    % Source depth
surface_z = 0;  % Surface location
H = -200;       % Depth
W = 10000;       % Maximum range
rho_air   = 1.225; % air density
rho_water = 1028;  % water density
rho_sand  = 1850;  % sand density

theta_aperture = 50;  % Aperture angle (in degree)
num_rays = 100;  % Number of rays to plot

source = [0, Z_src]; % Source location
absorption = 0.1;    % Absorption parameter (e.g., dB per meter)
% !!! Add frequency dependence !!!

% Updated sound speed profile
depths   = [   0, H];% -600, -4000]; %, -10,  -20,  -35,  -65, -100, -500 
c_values = [1500, 1300];%, 1520]; %, 1480, 1470, 1450, 1440, 1435, 1420
c = @(z) interp1(depths, c_values, z, 'linear'); % Interpolated sound speed function

figure()
plot(c_values, depths,'o', c(linspace(0,H,200)), linspace(0,H, 200))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Ray %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot rays and compute intensities
rays_r = cell(1, num_rays);
rays_z = cell(1, num_rays);
theta_values = linspace(-theta_aperture/2, theta_aperture/2, num_rays);

for i = 1:num_rays           % For each rays
    theta = theta_values(i) % Initial angle at the source

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%% Function %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = 0.01;      % Time step
    % Initialize the ray path
    r = [    0,         c(Z_src)*cosd(theta)*dt];
    z = [Z_src, Z_src + c(Z_src)*sind(theta)*dt];
    THETA = theta;
    while r(end) < W  % While the ray is in the rectangle
        
        % Actual conditions
        last_c    = c(z(end-1)); % Sound speed at last depth
        current_c = c(z(end  )); % Sound speed at current depth
        
        %g = (c(z(end)-1)-c(z(end)+1))/2;

        if abs((current_c/last_c)*cosd(theta)) < 1 % total transmission 
            theta = sign(theta)*acosd((current_c/last_c)*cosd(theta));
        else % total reflexion
            theta = -sign(theta)*acosd((last_c/current_c)*cosd(theta));
            H
        end

        THETA = [THETA, theta];

        % Calculate increments in r and z
        dr = current_c * cosd(theta) * dt;
        dz = current_c * sind(theta) * dt;
        
        % r = [r, r(end) + dr];
        % z = [z, z(end) + dz];
        
        % Check for reflection
        if z(end) + dz > 0          % on surface
            r = [r, r(end) - z(end)/tand(theta),   r(end) + dr ];
            z = [z,                            0, -(z(end) + dz)];
            surface_z
            theta=-theta;

        elseif z(end) + dz < H      % on bottom
            r = [r, r(end) + (H-z(end))/tand(theta),       r(end) + dr ];
            z = [z,                                H, -(z(end)+dz) + 2*H];
            surface_z-1
            theta=-theta;

        else    % Nominal case
            r = [r, r(end) + dr];
            z = [z, z(end) + dz];
        end

    end
    % figure()
    % plot(THETA)


    rays_r{i} = r;
    rays_z{i} = z;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Plot %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot setup
figure;
hold on;
xlabel('Range (m)');
ylabel('Depth (m)');
title('2D Ray Acoustic Modeling with Depth-Dependent Sound Speed');
axis([0 W H surface_z]);

% Plot surface and bottom
plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
plot([0, W], [H, H], 'k', 'LineWidth', 2);  % Bottom


for i = 1:num_rays
    r = rays_r{i};
    z = rays_z{i};
    plot(r, z)
end

hold off;
