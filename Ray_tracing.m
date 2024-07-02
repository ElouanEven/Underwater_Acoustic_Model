function [r, z, log_intensity, nbr_reflexion] = Ray_tracing(theta, SRC_z, H_d, c, W, absorption, dt)
    MAX_STEP = floor(W/(c(SRC_z)*dt) * 1.5); % Estimation of the max number of increments with theta_aperture < 90
    r = zeros(1,MAX_STEP);
    z = zeros(1,MAX_STEP);
    intensity = zeros(1,MAX_STEP);

    % Initialize the ray path and intensity
    r(1:2) = [    0,         c(SRC_z)*cosd(theta)*dt];
    z(1:2) = [SRC_z, SRC_z + c(SRC_z)*sind(theta)*dt];
    intensity(1:2) = [1, 10^(-absorption * c(z(end))*dt / 10)];
    count_ref = 0;
    i = 2;
    while (r(i) < W && i < MAX_STEP)  % While the ray is in the rectangle

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
    r = r(1:i);
    z = z(1:i);
    intensity = intensity(1:i);
    log_intensity = log10(intensity);
    nbr_reflexion = count_ref;
end