function [r, z, intensity, nbr_reflexion] = Ray_tracing(theta, SRC_z, I_0, H_d, c, W, absorption, dt, R_coef)
    % Ray_tracing - Compute one ray trace starting from the source with an angle theta under the influence of sound speed pressure with seabed and surface reflexion  
    %
    % Inputs:
    %    theta  - Angle from the ray at the source
    %    SRC_z  - Depth of the source
    %    I_0 - Ray intensity at the source
    %    H_d - Height of the seabed (depending on the range)
    %    c   - Sound speed profile (depending on the height)
    %    W   - Width max (max range)
    %    absorption  - Absorption (in dB/m)
    %    dt  - time step
    %    {c_sand, rho_sand, rho_water} - Environment condition (interesting for R_coef) !!! Not implemented yet !!!
    %    R_coef - Reflexion value for all grazing angle
    % 
    % Outputs:
    %    {r, z} - Ray coordonates
    %    intensity  - Ray intensity for each time step
    %    nbr_reflexion  - Number of reflexion (surface and bottom)
    %   
    %-----------------------------------------------------------------------------------------------
    delta = 0.01; % Width of the slope calcul for the seabed

    MAX_STEP = floor(W/(c(SRC_z)*dt) * 15); % Estimation of the max number of increments
    r = zeros(1,MAX_STEP);
    z = zeros(1,MAX_STEP);
    intensity = zeros(1,MAX_STEP);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize the ray path and intensity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dr_0 = c(SRC_z) * cosd(theta) * dt;
    dz_0 = c(SRC_z) * sind(theta) * dt;
    intensity_ref_0 = I_0 * 10^(-absorption * c(z(end))*dt / 10);

    r(1:2) = [    0,         dr_0];
    z(1:2) = [SRC_z, SRC_z + dz_0];
    intensity(1:2) = [I_0, intensity_ref_0]; % Intensity
    nbr_reflexion = 0;
    i = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Ray computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while (r(i) < W && i < MAX_STEP)  % While the ray is in the rectangle %%%%%%%%%%%%%%%%%%%%%%%%%%

        % Actual conditions
        last_c    = c(z(i-1)); % Sound speed at last depth
        current_c = c(z(i  )); % Sound speed at current depth

        if abs((current_c/last_c)*cosd(theta)) < 1 % total transmission     !!! Big problem on this condition !!!
            theta = sign(theta)*acosd((current_c/last_c)*cosd(theta));
        else % total reflexion
            theta = -sign(theta)*acosd((last_c/current_c)*cosd(theta));
        end
        
        % Calculate increments in r and z
        dr = current_c * cosd(theta) * dt;
        dz = current_c * sind(theta) * dt;
        
        % Check for reflection
        %% On surface  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if z(i) + dz > 0
            r(i+1:i+2) = [r(i) - z(i)/tand(theta), r(i) + dr ];
            z(i+1:i+2) = [0, -(z(i) + dz)];
            theta=-theta;
            intensity_ref = intensity(i) * 10^(-absorption * c(z(i))*dt / 10);
            intensity(i+1:i+2) = [intensity_ref, intensity_ref];
            nbr_reflexion = nbr_reflexion + 1;
            i=i+2;

        %% On bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif z(i) + dz < H_d(r(i) + dr)

            intensity_ref = intensity(i) * 10^(-absorption * c(z(i))*dt / 10);
            intensity_R = intensity_ref * 10^(-R_coef(theta) / 10);
            intensity(i+1:i+2) = [intensity_R, intensity_R];

            theta_bot = atand((H_d(r(i)+delta)-H_d(r(i)))/delta); % Bottom slop (calcul with 1m)
            r_bot = r(i) + (z(i)-H_d(r(i)))/(tand(abs(theta))+tand(theta_bot));     
            z_bot = H_d(r_bot);
            dist_ground = sqrt((r(i)+dr - r_bot)^2 + (z(i)+dz - z_bot)^2);
            theta = abs(theta)+2*theta_bot;

            if theta > 90 % No rays backward
                r = [r(1:i), r_bot];
                z = [z(1:i), z_bot];
                intensity = [intensity(1:i), intensity(i)];
                return
            end
            r_ref = r_bot + dist_ground*cosd(theta);
            z_ref = z_bot + dist_ground*sind(theta);
            r(i+1:i+2) = [r_bot, r_ref];
            z(i+1:i+2) = [z_bot, z_ref];

            nbr_reflexion = nbr_reflexion + 1;
            i=i+2;

        %% Nominal case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else    
            r(i+1) = r(i) + dr;
            z(i+1) = z(i) + dz;
            intensity(i+1) = intensity(i) * 10^(-absorption * c(z(i))*dt / 10);
            i=i+1;
        end
    end
    r = r(1:i);
    z = z(1:i);
    intensity = intensity(1:i);
end
