function [r, z] = Finding_Source(theta_detec, REC_r, REC_z_centered, H_d, c, W, dt)
    % Finding_Source - Ray tracing backward to find the source by using the correlation between time and space
    % 
    % Inputs:
    %    theta_detec - Angle detected on the receiver with the beamformer
    %    time_detec  - Depth of the source      !!! Not implemented yet !!!
    %    {REC_r, REC_z_centered} - Position of the center of the receiver
    %    H_d - Height of the seabed (depending on the range)
    %    c   - Sound speed profile (depending on the height)
    %    W   - Width max (max range)
    %    dt  - time step
    % 
    % Outputs:
    %    [r, z] - Ray coordonates
    %
    %-----------------------------------------------------------------------------------------------

    delta = 0.01; % Width of the slope calcul

    MAX_STEP = floor(W/(c(REC_z_centered)*dt) * 10); % Estimation of the max number of increments with theta_aperture < 90
    r = zeros(1,MAX_STEP);
    z = zeros(1,MAX_STEP);

    theta_detec = theta_detec + 180; % Rotation of 180 to reverse the rays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the ray path and intensity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dr_0 = c(REC_z_centered) * cosd(theta_detec) * dt;
    dz_0 = c(REC_z_centered) * sind(theta_detec) * dt;

    % Nominal condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r(1:2) = [REC_r,        REC_r + dr_0];
    z(1:2) = [REC_z_centered, REC_z_centered + dz_0];
    i = 2;
    theta = theta_detec;

    while (r(i) > -30 && i < MAX_STEP)  % While the ray is in the rectangle %%%%%%%%%%%%%%%%%%%%%%%%%%

        % Actual conditions
        last_c    = c(z(i-1)); % Sound speed at last depth
        current_c = c(z(i  )); % Sound speed at current depth
        
        if abs((current_c/last_c)*cosd(theta)) < 1 % total transmission     !!! Big problem on this condition !!!   %%%%%% Modif
            if theta < 180
                theta = acosd((current_c/last_c)*cosd(theta));
            else
                theta = 360 - acosd((current_c/last_c)*cosd(theta));
            end
        else % total reflexion
            if theta < 180
                theta = 360 - acosd((last_c/current_c)*cosd(theta));
            else
                theta = acosd((last_c/current_c)*cosd(theta));
            end
        end
        
        % Calculate increments in r and z
        dr = current_c * cosd(theta) * dt;
        dz = current_c * sind(theta) * dt;
        
        % Check for reflection
        % On surface  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if z(i) + dz > 0
            r(i+1:i+2) = [r(i) - z(i)/tand(theta), r(i) + dr ];
            z(i+1:i+2) = [0, -(z(i) + dz)];
            theta = 360-theta;
            i = i+2;

        % On bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif z(i) + dz < H_d(r(i) + dr)

            theta_bot = atand((H_d(r(i)-delta)-H_d(r(i)))/delta); % Bottom slop (calcul with 1m)
            r_bot = r(i) - (z(i)-H_d(r(i)))/(tand(abs(theta))+tand(theta_bot));     
            z_bot = H_d(r_bot);
            dist_ground = sqrt((r(i)+dr - r_bot)^2 + (z(i)+dz - z_bot)^2);
            theta = 360 - theta -2*theta_bot;

            r_ref = r_bot + dist_ground*cosd(theta);
            z_ref = z_bot + dist_ground*sind(theta);
            r(i+1:i+2) = [r_bot, r_ref];
            z(i+1:i+2) = [z_bot, z_ref];
            i=i+2;

        % Nominal case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else    
            r(i+1) = r(i) + dr;
            z(i+1) = z(i) + dz;
            i=i+1;
        end
    end
    r = r(1:i);
    z = z(1:i);
end
