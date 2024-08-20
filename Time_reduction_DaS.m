function [first_T_snap, last_T_snap] = Time_reduction_DaS(filtered_SIG, T_snap, nbr_rec, Nth_peak, Lfilt)

    PEAKS = [];
    LOCS = [];
    for i=1:nbr_rec
        [peaks, locs] = findpeaks(filtered_SIG(4*Lfilt:end,i), 'MinPeakDistance', 1000);
        PEAKS = [PEAKS, peaks'];
        LOCS = [LOCS, locs'];
    end
    
    [~, sortOrder] = sort(PEAKS, 'descend');
    
    % Apply the sorting order to all three vectors
    PEAKS = PEAKS(sortOrder);
    LOCS = LOCS(sortOrder);

    first_index = min(LOCS(1:Nth_peak));
    last_index = max(LOCS(1:Nth_peak));

    %%%%%%%%%%%%%%%% 
    
    first_T_snap = find(T_snap > first_index, 1) -3 
    last_T_snap = find(T_snap > last_index, 1) +3 +4
    if first_T_snap < 1
        first_T_snap = 1;
    end
    if last_T_snap > length(T_snap)
        last_T_snap = length(T_snap)-1;
    end
    
    % Plot time signal
    figure;
    for i=1:nbr_rec
        plot(filtered_SIG(:,i));
        hold on;
    end
    % Plot all peaks
    plot(LOCS(1:Nth_peak)+4*Lfilt, PEAKS(1:Nth_peak), 'rv', 'MarkerFaceColor', 'r');  % Red inverted triangle for peaks
    hold off;
end
