function [ID_ray, delay_t, delay_i] = Cross_receiver(rays_r, rays_z, rays_log_intensity, REC_r, REC_z, width, num_rays, dt, nbr_reflexion)
    ID_ray = [];
    delay_t = [];
    delay_i = [];
    for i = 1:num_rays
        for j = 1:length(rays_r{1,i})-1
            if (rays_r{1,i}(j) < REC_r) && (rays_r{1,i}(j+1) > REC_r)
                A = (rays_z{1,i}(j+1)-rays_z{1,i}(j))/(rays_r{1,i}(j+1)-rays_r{1,i}(j));
                B = rays_z{1,i}(j) - A*rays_r{1,i}(j);
                f_REC_r = A*REC_r + B;
                if (f_REC_r > REC_z) && (f_REC_r < REC_z+width)
                    ID_ray = [ID_ray, i];
                    delay_t = [delay_t, dt*(j - nbr_reflexion(i))];
                    delay_i = [delay_i, rays_log_intensity{1,i}(j)];
                end
            end
        end
    end
end