function [] = Plot_All_Rays(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, ZONE_width)
    REC_z_centered = mean(REC_z); % Center of the receiver
    figure;
    hold on;
    
    % Adjusting font size for labels and title
    xlabel('Range (m)', 'FontSize', 30, 'FontWeight', 'bold');
    ylabel('Depth (m)', 'FontSize', 30, 'FontWeight', 'bold');
    title('Acoustic rays emitted by the source', 'FontSize', 30, 'FontWeight', 'bold');
    
    % Plotting the surface and bottom
    axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
    plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
    plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'r', 'LineWidth', 2);  % Bottom
    
    % Plotting all rays without color
    for i = 1:length(rays_r)
        r = rays_r{i};
        z = rays_z{i};
        line(r, z); % Plot all rays
    end
    
    % Plotting the zone and receiver
    plot([REC_r; REC_r], [REC_z_centered-ZONE_width/2; REC_z_centered+ZONE_width/2], 'b', 'LineWidth', 8);  % Zone
    plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
    
    hold off;
    % Set background color to white
    set(gcf, 'Color', 'w');

    % Adjusting general font size and tick size for the axes
    set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLength', [0.02, 0.02]);
    
    % Optional: Customize tick direction and other properties if needed
    set(gca, 'TickDir', 'out');  % Tick direction can be 'in' or 'out'

end
