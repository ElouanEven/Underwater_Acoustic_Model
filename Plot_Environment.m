function [] = Plot_Environment(H_values, H_d, distance, W, surface_z, SRC_z, REC_r, REC_z, width)
    figure;
    hold on;
    axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
    plot(distance, H_values, 'o', 'LineWidth', 2);
    plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'b', 'LineWidth', 2);
    plot(0, SRC_z, 'o', 'LineWidth', 5);                % Source
    plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
    plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
    % Plot surface and seabed
    xlabel('Distance from source (m)', 'LineWidth', 2);
    ylabel('Altitude (m)', 'LineWidth', 2);
    title('Bathymetrie', 'LineWidth', 2);
    legend('Data', 'Interpolation', 'Source', 'Receiver', 'LineWidth', 2);
    hold off;
end