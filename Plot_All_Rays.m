function [] = Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z)
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
    
    
    % For all rays without color
    for i = 1:length(rays_r)
        r = rays_r{i};
        z = rays_z{i};
        line(r, z); % Plot all rays
    end
    
    hold off;
end