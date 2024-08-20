function [max_t, max_phi] = Find_local_maxima(ENERGY, nbr_maxima, time_width)

    localMaxima = imregionalmax(ENERGY); % Find local maxima
    [max_t_tot, max_phi_tot] = find(localMaxima); % Extract the coordinates of the local maxima
    max_values = ENERGY(sub2ind(size(ENERGY), max_t_tot, max_phi_tot)); % Get the values of the local maxima
    [~, sortedIndices] = sort(max_values, 'descend'); % Sort the maxima by their values in descending order
    
    top_Indices = sortedIndices(1:nbr_maxima); % Select the top n maxima
    max_t = max_t_tot(top_Indices);
    max_phi = max_phi_tot(top_Indices);

    % Tests if max unique on a time
    indexToRemove = [];
    for i =1:nbr_maxima-1
        for j=i+1:nbr_maxima
            if abs(max_t(i)-max_t(j)) < time_width
                indexToRemove = [indexToRemove, j];
            end
        end
    end
    indexToRemove = unique(indexToRemove);
    % Remove the value if found
    if ~isempty(indexToRemove)
        max_t(indexToRemove) = [];
        max_phi(indexToRemove) = [];
    end
end