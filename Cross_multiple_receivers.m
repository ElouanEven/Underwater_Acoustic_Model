function [ID_ray, delay_t, delay_i] = Cross_multiple_receivers(rays_r, rays_z, rays_intensity, REC_r, REC_z, width, dt, nbr_reflexion, c_REC_z, absorption)
    % Cross_multiple_receivers - Verifies if the ray is crossing one of the receiver 
    %
    % Inputs:
    %    {rays_r, rays_z} - Coordonate of each point of the rays
    %    rays_intensity - Ray intensity for each time step
    %    {REC_r, REC_z} - Coordonate of the receiver
    %    width          - Receiver width
    %    dt             - Ray intensity at the source
    %    nbr_reflexion  - Number of reflexion (surface and bottom)
    %    c_REC_z        - Sound speed at the source (fixed)
    %    absorption     - Absorption (in dB/m)
    % 
    % Outputs:
    %    ID_ray - # of the ray
    %    {delay_t, delay_i} - delay and intensity of the ray on the receiver
    %
    %-----------------------------------------------------------------------------------------------

    nbr_rec = length(REC_z);
    ID_ray = cell(1, nbr_rec);
    delay_t = cell(1, nbr_rec);
    delay_i = cell(1, nbr_rec);
    for i = 1:length(rays_r) % For each rays
        for j = 1:length(rays_r{1,i})-1 % For each point of this ray
            if (rays_r{1,i}(j) < REC_r) && (rays_r{1,i}(j+1) > REC_r) % Select the two points arround the range of the receiver
                % Finding the equation of the line A*X+B
                A = (rays_z{1,i}(j+1)-rays_z{1,i}(j))/(rays_r{1,i}(j+1)-rays_r{1,i}(j));
                B = rays_z{1,i}(j) - A*rays_r{1,i}(j);
                f_REC_r = A*REC_r + B;
                for k =1:nbr_rec % For each receiver
                    if (f_REC_r > REC_z(k)) && (f_REC_r < REC_z(k)+width) % If the line between the two points cross the k-th receiver
                        
                        dist_ray_rec = sqrt((REC_r-rays_r{1,i}(j))^2 + (f_REC_r-rays_z{1,i}(j))^2); % Distance between the last point of the ray and the k-th receiver
                        
                        ID_ray {1,k} = [ ID_ray{1,k}, i];
                        delay_t{1,k} = [delay_t{1,k}, dt*(j - nbr_reflexion(i)) + dist_ray_rec/c_REC_z]; 
                        delay_i{1,k} = [delay_i{1,k}, rays_intensity{1,i}(j)*10^(-absorption * dist_ray_rec / 10)];
                    end
                end
            end
        end
    end
end
