function [] = Plot_Receiver_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, ID_ray, delay_i)
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
    % Color mapping based on log intensity
    colormap('jet');
    JET = jet(256);
    
    if length(ID_ray) == 1 % If only one ray touch the receiver
        line(rays_r{1}, rays_z{1});  % Apply color to line
    else
        j=0;
        for i = ID_ray
            r = rays_r{i};
            z = rays_z{i};
            j=j+1;
            color_idx = 256 - floor((delay_i(j) - min(delay_i)) / (max(delay_i) - min(delay_i)) * 255);
            color = JET(color_idx, :);  % Get color from colormap based on log intensity
            line(r, z, 'Color', color);  % Apply color to line
            
        end
        % Add color bar for intensity
        caxis([min(delay_i), max(delay_i)]);
        colorbar;
        ylabel(colorbar, 'Log Intensity');
    end
    hold off;
end