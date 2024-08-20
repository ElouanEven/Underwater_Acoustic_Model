function [r_inter, z_inter, t_match, Best_time, R_SRC, Z_SRC] = Find_Intersection(rays_r_det, rays_z_det, dt, list_time_detected, list_angle_detected)
    % Trouver toutes les intersections
    r_inter = [];
    z_inter = [];
    t_match = [];
    Best_time = 10;
    R_SRC = 0;
    Z_SRC = 0;
    % For each couple of rays
    for i = 1:length(rays_r_det)-1
        for j = i+1:length(rays_r_det)
            r1 = rays_r_det{1,i};
            z1 = rays_z_det{1,i};
            r2 = rays_r_det{1,j};
            z2 = rays_z_det{1,j};

            n1 = length(r1);
            n2 = length(r2);
            % For each couple of segment of both rays
            for k = 100:n1-1
                for l = 1:n2-1
                    if abs(r1(k)-r2(l))<0.5 && abs(z1(k)-z2(l))<0.5% && abs(list_angle_detected(i)-list_angle_detected(j)) > 3 % Test with a minimum of angle between two comparison
                        [ri, zi] = polyxpoly([r1(k) r1(k+1)], [z1(k) z1(k+1)], ...
                                             [r2(l) r2(l+1)], [z2(l) z2(l+1)]);
                        if ~isempty(ri)
                            r_inter = [r_inter; ri];
                            z_inter = [z_inter; zi];
                            delta_t = abs(k*dt - l*dt);
                            Time_expected = abs(list_time_detected(i) - list_time_detected(j));
                            t_match = [t_match, abs(Time_expected-delta_t)];
                            if abs(Time_expected-delta_t)<Best_time
                                Best_time = abs(Time_expected-delta_t);
                                R_SRC = ri;
                                Z_SRC = zi;
                            end
                        end
                    end
                end
            end
        end
    end
    if R_SRC==0 && Z_SRC==0
        disp("!!!!! No intersection between rays !!!!!");
    end
    [~, sortOrder] = sort(t_match);
    
    % Apply the sorting order to all three vectors
    r_inter = r_inter(sortOrder);
    z_inter = z_inter(sortOrder);
    t_match = t_match(sortOrder);
end