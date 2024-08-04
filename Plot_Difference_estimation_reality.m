function [] = Plot_Difference_estimation_reality(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, ID_ray, rays_r_det, rays_z_det)
    nbr_rec = length(REC_z);
    figure;
    hold on;
    xlabel('Range (m)');
    ylabel('Depth (m)');
    title('Rac SSP receiver');
    axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
    % Plot surface and bottom
    plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
    plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'r', 'LineWidth', 2);  % Bottom
    plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
    % Color mapping based on intensity
    colormap('jet');
    JET = jet(256);
    for k = 1:nbr_rec
        color = floor((k / nbr_rec) * 255) + 1;
        for i = ID_ray{1,k}
            line(rays_r{1,i}, rays_z{1,i}, 'Color', JET(color, :));
        end
    end
    
    % For all rays detected by the receiver
    for i = 1:length(rays_r_det)
        r = rays_r_det{i};
        z = rays_z_det{i};
        line(r, z, 'LineWidth', 3, 'LineStyle', ':', 'Color', 'r'); % Plot all rays
    end
    hold off;
end