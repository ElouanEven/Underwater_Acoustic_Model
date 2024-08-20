function [ENERGY, list_angle_detected, list_time_detected] = Beamformer_DaS(filtered_SIG, REC_z, Lfilt, F, PHI, c_iso, T_0, T_snap)
    % Beamformer_DaS - Applying the Delay and Sum Beamformer to find the direction of arrival over time to the input signal
    %
    % Inputs:
    %    filtered_SIG - Signal arriving on each sensor
    %    REC_z     - Position of each sensor
    %    Lfilt     - Frequency resolution
    %    F         - Frequency (list)
    %    PHI       - Angle (list)
    %    c_iso     - Sound speed supposed homogeneous on the profile !!!Wrong Hypothesis!!!
    %    T_0       - Time vector
    %    T_snap    - Time for each start snapshot with overlap
    %    threshold - Min amplitude to consider a direction
    % 
    % Outputs:
    %    ENERGY - Energy of the directivity pattern
    %    list_angle_detected - Angle of arrival on the receiver (multiple with time)
    %    list_time_detected - Time associated to direction of arrival
    %
    %-----------------------------------------------------------------------------------------------

    N_phi = length(PHI);
    nbr_rec = length(filtered_SIG(1,:));
    
    HannW = meshgrid(hann(Lfilt), zeros(nbr_rec,1)); % Hanning window to multiply to each x_snap sensor
    
    list_angle_detected = [];
    list_time_detected = [];

    Nth_peak = 30;

    nbr_maxima = 7;
    time_width = 2;

    [first_T_snap, last_T_snap] = Time_reduction_DaS(filtered_SIG, T_snap, nbr_rec, Nth_peak, Lfilt); % Reduction of the time of interest


    ENERGY = zeros(last_T_snap-first_T_snap,N_phi);

    for t_snap = first_T_snap:last_T_snap %1:length(T_snap)
        SIG_snap = filtered_SIG(T_snap(t_snap)+1:T_snap(t_snap)+Lfilt,:);
        SIG_snap_window = HannW' .* SIG_snap; % Application of the hanning window on the snapshot signal
        FFT_x = fft(SIG_snap_window);         % Fourier transform
        Y = zeros(length(F), length(PHI));    % Initializing the output array
        
        for F_i = 1:length(F)
            for PHI_i = 1:length(PHI)
                w = zeros(nbr_rec,1); % Initializing weights for all sensors 
                for rec_i = 1:nbr_rec
                    u0 = [cosd(PHI(PHI_i)); 0; 0];          % Direction normalised
                    w(rec_i) = 1/nbr_rec * exp(-1i*2*pi/c_iso.*F(F_i)*[REC_z(rec_i)-REC_z(1) 0 0]*u0); % Weight equation
                end
                
                Y(F_i, PHI_i) = FFT_x(F_i,:) * w; % Output signal
            end
        end
        
        Y_abs = abs(Y);     % Absolute signal
        Energy = Y_abs.^2;  % Energy
    
        ENERGY(t_snap-first_T_snap+1,:) = mean(Energy,1);
    end

    
    [max_t, max_phi] = Find_local_maxima(ENERGY, nbr_maxima, time_width);

    list_time_detected = T_0(round(T_snap(max_t)));
    list_angle_detected = PHI(max_phi)-90;
end
