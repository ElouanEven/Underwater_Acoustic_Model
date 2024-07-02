function [] = Plot_All_Rays_Intensity(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, rays_log_intensity, log_intensity_min, log_intensity_max)
    figure;
    hold on;
    xlabel('Range (m)');
    ylabel('Depth (m)');
    title('Rac SSP all rays');
    axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
    % Plot surface and bottom
    plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
    plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'r', 'LineWidth', 2);  % Bottom
    plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
    % Color mapping based on log intensity
    colormap('jet');
    JET = jet(256);
    
    % For all rays with color intensity
    for i = 1:num_rays
        r = rays_r{1,i};
        z = rays_z{1,i};
        log_intensity = rays_log_intensity{1,i};

       for j = 1:length(r)-1
            color_idx = floor((log_intensity(j) - log_intensity_min) / (log_intensity_max - log_intensity_min) * 255) + 1;
            color = JET(color_idx, :);  % Get color from colormap based on log intensity
            line([r(j), r(j+1)], [z(j), z(j+1)], 'Color', color);  % Apply color to line
        end
    end
    
    
    % Add color bar for intensity
    caxis([log_intensity_min, log_intensity_max]);
    colorbar;
    ylabel(colorbar, 'Log Intensity');
    
    hold off;
end