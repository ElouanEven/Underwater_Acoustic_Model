function [ID_ray_zone] = Cross_zone(rays_r, rays_z, REC_r, REC_z_centered, ZONE_width)
    % Cross_zone - Check which rays are crossing the zone
    %
    % Inputs:
    %    {rays_r, rays_z} - Ray coordonates
    %    {REC_r, REC_z_centered} - Position of the center of the zone
    %    ZONE_width - Width of the zone
    % 
    % Outputs:
    %    ID_ray_zone - List of # of rays crossing the zone
    %   
    %-----------------------------------------------------------------------------------------------
    ID_ray_zone = [];
    for i = 1:length(rays_r) % For each ray
        for j = 1:length(rays_r{1,i})-1 % For each segment of the ray
            if (rays_r{1,i}(j) < REC_r) && (rays_r{1,i}(j+1) > REC_r) % Select the two points arround the range of the zone
                % Finding the equation of the segment A*X+B
                A = (rays_z{1,i}(j+1)-rays_z{1,i}(j))/(rays_r{1,i}(j+1)-rays_r{1,i}(j));
                B = rays_z{1,i}(j) - A*rays_r{1,i}(j);
                f_REC_r = A*REC_r + B;
                if (f_REC_r > REC_z_centered-ZONE_width/2) && (f_REC_r < REC_z_centered+ZONE_width/2) % Test if the z value cross the zone
                    ID_ray_zone = [ID_ray_zone, i];
                end
            end
        end
    end
end