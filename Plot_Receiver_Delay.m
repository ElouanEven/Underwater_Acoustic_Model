function [] = Plot_Receiver_Delay(nbr_reflexion, ID_ray, delay_t, delay_i)
    JET = jet(256);
    
    if max(nbr_reflexion(ID_ray)) == min(nbr_reflexion(ID_ray))
        Delay_color = ones(1,length(delay_t));
    else
        Delay_color = floor((nbr_reflexion(ID_ray) - min(nbr_reflexion(ID_ray))) / (max(nbr_reflexion(ID_ray)) - min(nbr_reflexion(ID_ray))) * 255) + 1;
    end
    
    figure;
    hold on;
    for i = 1:length(delay_t)
        plot(delay_t(i), delay_i(i), 'o', 'Color', JET(Delay_color(i), :), 'MarkerSize', 8, 'LineWidth', 2);
    end
    hold off;
end