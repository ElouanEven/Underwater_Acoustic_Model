function [micIRs, micTFs] = ArrayResponse(U_doa, pn, U_orient, fDir_handle, Nfft, fs)
% For now we will consider only one source so one U_doa vector
% We won't consider the orientation of the sensors

%% %%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%

Nmics = size(pn,1); % Number of microphones
c = 340;           % speed of sound
K = Nfft/2+1;       % Number of frequency
f = (0:K)*fs/Nfft;  % frequency vector

%%%%%%%%%% COMPUTATIONS %%%%%%%%%%
% orientation of the microphones is assumed to be radial from the origin
%U_orient = pn./sqrt(sum(pn.^2,2));

% unit vectors pointing to the evaluation points
U_eval = meshgrid(U_doa,ones(1,Nmics));

dcosAngleU = dot(U_eval, pn, 3);

% create TFs for each microphone (delay)
micTFs = zeros(K, Nmics);
for kk = 1:K
    omega = 2*pi*f(kk);
    tempTF = exp(1i*(omega/c)*dcosAngleU(:,3)); % !!! Only angle z !!!
    micTFs(kk,:) = tempTF.';
end
size(micTFs)

% create IRs for each microphone

tempTF = micTFs(:,:);
tempTF(end,:);
tempTF(end,:) = abs(tempTF(end,:));
tempTF = [tempTF; conj(tempTF(end-1:-1:2,:))];
micIRs(:,:) = ifft(tempTF);
micIRs(:,:) = fftshift(micIRs(:,:), 1);

