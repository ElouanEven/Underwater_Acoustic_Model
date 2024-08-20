function [] = Plot_Difference_estimation_reality(H_d, W, surface_z, REC_r, REC_z, width, rays_r, rays_z, ID_ray, rays_r_det, rays_z_det, r_inter, z_inter)
    nbr_rec = length(REC_z);
    x = linspace(0, W, 200);
    y = H_d(x);
    maxY = min(H_d(linspace(0, W, 200))) - 2;
    brown = [0.65, 0.16, 0.16];
    colormap('jet');
    JET = jet(256);

    figure;
    hold on;
    
    % Adjusting font size for labels and title
    xlabel('Range (m)', 'FontSize', 30, 'FontWeight', 'bold');
    ylabel('Depth (m)', 'FontSize', 30, 'FontWeight', 'bold');
    title('Acoustic ray propagation and estimation of the source', 'FontSize', 30, 'FontWeight', 'bold');
    
    % Plotting the surface and bottom
    axis([0 W min(H_d(linspace(0, W, 200))) surface_z]);
    plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
    plot(linspace(0, W, 200), H_d(linspace(0, W, 200)), 'Color', brown, 'LineWidth', 2);  % Bottom
    fill([x fliplr(x)], [y repmat(maxY, 1, length(x))], brown, 'FaceAlpha', 0.5, 'EdgeColor', 'none');  % Filled area for bottom
    plot([REC_r; REC_r], [REC_z; REC_z+width], 'r', 'LineWidth', 4);  % Receiver
    
    % Plotting all straight rays with the same color and increased line width
    ray_line = line(nan, nan, 'Color', 'g', 'LineWidth', 1);  % Placeholder for legend
    for k = 1:nbr_rec
        color = floor((k / nbr_rec) * 255) + 1;
        for i = ID_ray{1,k}
            line(rays_r{1,i}, rays_z{1,i}, 'Color', JET(color, :));
        end
    end
    
    % Plotting all rays detected by the receiver with dotted lines
    ray_est_line = line(nan, nan, 'LineStyle', ':', 'LineWidth', 3, 'Color', 'r');  % Placeholder for legend
    for i = 1:length(rays_r_det)
        r = rays_r_det{i};
        z = rays_z_det{i};
        r_l=line(r, z, 'LineWidth', 3, 'LineStyle', ':', 'Color', 'r');  % Rays estimation
    end
    
    % Plotting possible source (cross) and creating placeholder for the legend
    source_plot = plot(r_inter, z_inter, 'o', 'LineWidth', 5, 'Color', 'b');  % Possible source
    
    hold off;
    
    % Creating the legend with all necessary elements
    legend([ray_line, ray_est_line, source_plot], 'Rays from source', 'Rays estimation', 'Possible source', 'FontSize', 30, 'Location', 'northeast');
    
    % Adjusting general font size and tick size for the axes
    set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLength', [0.02, 0.02]);
    
    % Set background color to white
    set(gcf, 'Color', 'w');
    
    % Optional: Customize tick direction and other properties if needed
    set(gca, 'TickDir', 'out');  % Tick direction can be 'in' or 'out'


end
