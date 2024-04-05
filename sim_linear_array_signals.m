% For now we will consider only one source
clear
close all
clc

%% %%%%%%%%%%%%%%%%%%---- Array geometry -----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% Uniform linear array
Q = 20;        % number of sensors
q = (0:Q-1)'; % index
d = 0.5;      % inter-sensor distance

zn = ((q-(Q-1)/2)*d); % sensor positions on x-axis
xn = zeros(Q,1);
yn = zeros(Q,1);
pn = [xn yn zn];      % microphone positions

% source direction
doa_elev_deg = 20;
doa_elev = doa_elev_deg*pi/180; % elevation follows Matlab definition [-90, 90]
doa_azi_deg = 0;
doa_azi= doa_azi_deg*pi/180; 

[x,y,z] = sph2cart(doa_azi, doa_elev, 1); % cartesian coordinates
U_doa = [x y z]; % unit vector

% Impulse response parameters
Fs = 2000;      % Sampling frequency
Lfilt = 800;   % Max Frequency

% simulate array response using getArrayResponse()
fDirectivity = @(angle) 1; % response of omnidirectional sensor
[h_mic, H_mic] = getArrayResponse(U_doa, pn, [], fDirectivity, Lfilt, Fs); % sensor orientation irrelevant in this case

fprintf('Array geometry: %.4f seconds\n', toc);

%% %%%%%%%%%%%%%%%%%%---- Static source -----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

dur = 1; % duration in seconds
N = round(Fs*dur); % duration in samples

Ps = 5; % signal power 
sig = sqrt(Ps)*randn(N,1);


xs = zeros(N,Q); % array signals
for i=1:Q
    xs(:,i) = fftfilt(h_mic(:,i),sig);
end

fprintf('Static source: %.4f seconds\n', toc);

%% %%%%%%%%%%%%%%%%%---- White noise ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
SNR = 0;

Pn = Ps/(10^(SNR/10)); % noise power at each sensor

sn = sqrt(Pn)*randn(N,Q);

x_tot = xs + sn;

fprintf('White noise: %.4f seconds\n', toc);
%% %%%%%%%%%%%%%%%%%---- Beamforming ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

T = 0:1/Fs:dur-1/Fs;
F = (0:Lfilt/2-1)*Fs/Lfilt;
N_phi = 360;
PHI = 0:2*pi/N_phi:2*pi-2*pi/N_phi;
c = 340;

[F_mesh, PHI_mesh, Q_mesh] = meshgrid(F, PHI, q);

W = 1/Q * exp(-1j*2*pi/c*d*cos(PHI_mesh).*F_mesh.*Q_mesh);

size_W = size(W);
size_x = size(x_tot);
W_reshaped = reshape(W, size_W(1)*size_W(2), size_W(3)); % reshape W to be [size(F*PHI)xQ]

% Perform matrix multiplication
Y_reshaped = W_reshaped * x_tot'; % dimensions will be [size(F*PHI)xsize(x_tot)]

% Reshape the result back to the desired dimensions
Y = reshape(Y_reshaped, size_W(1), size_W(2), size_x(1));

Y_abs = abs(Y);
Energy = Y.^2;

fprintf('Beamforming: %.4f seconds\n', toc);

%% %%%%%%%%%%%%%%%%%---- Visualisation ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHI_mesh_deg = PHI_mesh*180/pi;

[T_mesh,PHIT_mesh_deg] = meshgrid(T, PHI*180/pi);

% X_plot = [0 0 T(end) T(end)]; % X coordinates
% Y_plot = [doa_elev_deg doa_elev_deg doa_elev_deg doa_elev_deg]; % Y coordinates
% Z_plot = [min(Y_abs(:,1,:)) max(Y_abs(:,1,:)) min(Y_abs(:,1,:)) max(Y_abs(:,1,:))]; % Z coordinates 
% [X_plot2, Y_plot2] = meshgrid(X_plot, Y_plot);

% Array geometry
figure (1)
plot3(pn(:,1),pn(:,2),pn(:,3),'ko','MarkerSize',8)
grid on
title('Position of sensors')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
set(gca, 'fontsize', 20, 'linewidth',2)
set(findall(gcf, 'Type', 'Line'),'LineWidth',2); 
view(120,30)

figure(2)
plot(T, x_tot);
title('Signal received from the sensors');
xlabel('Time (s)');
ylabel('Amplitude');
set(gca, 'fontsize', 20, 'linewidth',2)
grid on;

figure(3)
mesh(F_mesh(:,:,1), PHI_mesh_deg(:,:,1), abs(Energy(:,:,1)));%  Y_abs(:,:,1));
colormap(jet);
xlabel('Frequency (Hz)');
ylabel('Phi (deg)');
zlabel('Amplitude');
title('3D Plot with Colormap of Y_{tot}');
colorbar;

figure(4)
mesh(T_mesh(:,:,1), PHIT_mesh_deg(:,:,1), squeeze(Y_abs(:,1,:)));
colormap(jet);
xlabel('Time (s)');
ylabel('Phi (deg)');
zlabel('Amplitude');
title('3D Plot with Colormap of Y_{tot}');
colorbar;




