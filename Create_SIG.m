function [SIG] = Create_SIG(nbr_rec, len_0, ID_ray, sig, Fs, delay_t, delay_i, rho_water, c_water_mean)
    % Create_SIG - Create the signal on each receiver
    %
    % Inputs:
    %    nbr_rec        - Number of receiver on the receiver array
    %    len_0          - Approximation of the total time length vector
    %    ID_ray         - # of the ray
    %    sig            - Signal from the source
    %    Fs             - Sampling frequency (Hz)
    %    {delay_t, delay_i} - delay and intensity of the ray on the receiver
    %    rho_water      - Water density
    %    c_water_mean   - Approximation of sound speed in water
    % 
    % Outputs:
    %    SIG - Signal on each receiver
    %
    %-----------------------------------------------------------------------------------------------

    SIG = zeros(len_0, nbr_rec);
    for k = 1:nbr_rec % For all receivers
        for i = 1:length(ID_ray{1,k}) % For all rays arriving on the receiver k
            i_timedelay = round(delay_t{1,k}(i)*Fs);
            i_start = i_timedelay;
            i_end = i_timedelay+(length(sig));
            sig_delay = sig*sqrt(delay_i{1,k}(i) * 2*rho_water*c_water_mean); % Do we have to use the c_water_mean or the c(rec_z) ???
            if i_end < len_0
                SIG(i_start:i_end-1,k) = SIG(i_start:i_end-1,k) + sig_delay';
            else
                fprintf("Issue with ray : %.5f \n", delay_t{1,k}(i)); % The ray didn't reach the receiver before the last step ==> Increase the value of MAX_STEP
            end
        end
    end

end