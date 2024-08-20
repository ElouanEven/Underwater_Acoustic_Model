function [theta_values_zone] = Theta_Improved(ID_ray_zone, theta_values, theta_resolution, coef_improvment_Dthteta)
    interesting_angles = theta_values(ID_ray_zone); % All angles from the source leading to cross the zone
    theta_res_vec = linspace(-theta_resolution, theta_resolution, coef_improvment_Dthteta*2+1); % Improving the resolution angle
    theta_values_zone = [];
    for i = 1:length(interesting_angles)
        theta_values_zone = [theta_values_zone, interesting_angles(i) + theta_res_vec];
    end
    theta_values_zone = unique(theta_values_zone);
end