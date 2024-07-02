function [] = Plot_SSP(c_values, c, depths, H, surface_z)
    figure;
    hold on;
    plot(c_values, depths, 'o', 'LineWidth', 2);
    plot(c(linspace(0, H, 200)), linspace(0, H, 200), 'LineWidth', 2); 
    % Plot surface and bottom
    plot([min(c(linspace(0, H, 200))), max(c(linspace(0, H, 200)))], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
    plot([min(c(linspace(0, H, 200))), max(c(linspace(0, H, 200)))], [H, H], 'k', 'LineWidth', 2);  % Bottom
    xlabel('Sound speed (m/s)', 'LineWidth', 2);
    ylabel('Depth (m)', 'LineWidth', 2);
    title('Sound Speed Profile', 'LineWidth', 2);
    legend('Data', 'Interpolation', 'LineWidth', 2);
    hold off;
end