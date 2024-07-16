function [ID_ray, delay_t, delay_i] = Cross_multiple_receivers(rays_r, rays_z, rays_intensity, REC_r, REC_z, width, dt, nbr_reflexion, c_REC_z, absorption)
    nbr_rec = length(REC_z);
    ID_ray = cell(1, nbr_rec);
    delay_t = cell(1, nbr_rec);
    delay_i = cell(1, nbr_rec);
    theta = [];
    for i = 1:length(rays_r)
        for j = 1:length(rays_r{1,i})-1
            if (rays_r{1,i}(j) < REC_r) && (rays_r{1,i}(j+1) > REC_r)
                % Finding the equation of the line A*X+B
                A = (rays_z{1,i}(j+1)-rays_z{1,i}(j))/(rays_r{1,i}(j+1)-rays_r{1,i}(j));
                B = rays_z{1,i}(j) - A*rays_r{1,i}(j);
                f_REC_r = A*REC_r + B;
                for k =1:nbr_rec
                    if (f_REC_r > REC_z(k)) && (f_REC_r < REC_z(k)+width) 
                        
                        dist_ray_rec = sqrt((REC_r-rays_r{1,i}(j))^2 + (f_REC_r-rays_z{1,i}(j))^2); % Distance between the  last ray and the receiver

                        ID_ray{1,k} = [ID_ray{1,k}, i];
                        delay_t{1,k} = [delay_t{1,k}, dt*(j - nbr_reflexion(i)) + dist_ray_rec/c_REC_z];
                        delay_i{1,k} = [delay_i{1,k}, rays_intensity{1,i}(j)*10^(-absorption * dist_ray_rec / 10)];
                        theta = [theta, atand(A)];
                    end
                end
            end
        end
    end
    edges = -90:1:90;
    figure;
    histogram(theta, 'BinEdges', edges);
    xlabel('Data Value');
    ylabel('Frequency');
    title('Histogram with Specified Bin Edges and Number of Bins');
end
