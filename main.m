clc; clearvars; close all;

%% -------------------------------------------------------------------------
% 1. Setup system and model parameters.
% -------------------------------------------------------------------------
% System parameters
EIRP_dBm = 60; % transmit array EIRP
P_noise_dBm = -68; % receive array integrated noise power (includes amplification)

% Model parameters: location and scale parameters for mu
G_dB = -129.00; % location parameter
xi = 0.502; % scale parameter

% Model parameters: cluster centers and angular spreads for H_bar
aod_az_list = [-174 126 -118 126]; % AoD azimuth
aod_el_list = [0 0 0 0]; % AoD elevation
aoa_az_list = [-122 -122 -122 118]; % AoA azimuth
aoa_el_list = [0 0 0 0]; % AoA elevation
spread_az = 4; % angular spread in azimuth
spread_el = 3; % angular spread in elevation
spread_az_res = 1; % angular spread resolution in azimuth
spread_el_res = 1; % angular spread resolution in elevation

% Model parameters: estimator and variance parameters for sigma^2
alpha = -0.733; % slope
beta = 42.53; % bias
nu_squared = 126.091; % variance

%% -------------------------------------------------------------------------
% 2. Create array objects.
% -------------------------------------------------------------------------
% planar array size
M = 16; % number of rows
N = 16; % number of columns

% create transmit and receive arrays (custom array object)
atx = array.create(M,N); % transmit array
arx = array.create(M,N); % receive array

% number of transmit and receive antennas
Nt = M * N;
Nr = M * N;

%% -------------------------------------------------------------------------
% 3. Define transmit and receive codebook steering directions.
% -------------------------------------------------------------------------
% transmit directions (e.g., transmit codebook)
tx_dir_az = flip([-56:8:56]); % flip is not necessary
tx_dir_el = flip([-8:8:8]);

% each transmit az-el pair
tx_dir_az_el_deg = [repelem(tx_dir_az(:),length(tx_dir_el)),repmat(tx_dir_el(:),length(tx_dir_az),1)];
num_tx = length(tx_dir_az_el_deg(:,1));

% receive directions (e.g., receive codebook)
rx_dir_az = tx_dir_az; % assume same for simplicity
rx_dir_el = tx_dir_el;

% each receive az-el pair
rx_dir_az_el_deg = [repelem(rx_dir_az(:),length(rx_dir_el)),repmat(rx_dir_el(:),length(rx_dir_az),1)];
num_rx = length(rx_dir_az_el_deg(:,1));

%% -------------------------------------------------------------------------
% 4. Transmit and receive beamforming codebooks.
% -------------------------------------------------------------------------
% array response vectors in each steering direction
Atx = atx.get_array_response(tx_dir_az_el_deg(:,1)*pi/180,tx_dir_az_el_deg(:,2)*pi/180);
Arx = arx.get_array_response(rx_dir_az_el_deg(:,1)*pi/180,rx_dir_az_el_deg(:,2)*pi/180);

% transmit and receive codebooks use conjugate beamforming
F = Atx;
W = Arx;


% ensure beams in F are normalized
for idx_tx = 1:num_tx
    f = F(:,idx_tx);
    F(:,idx_tx) = f ./ norm(f,2) .* sqrt(Nt);
end

% ensure beams in W are normalized
for idx_rx = 1:num_rx
    w = W(:,idx_rx);
    W(:,idx_rx) = w ./ norm(w,2) .* sqrt(Nr);
end

%% -------------------------------------------------------------------------
% 5. Construct channel matrix (H bar).
% -------------------------------------------------------------------------
H = zeros(Nr,Nt);
for idx_clust = 1:length(aod_az_list)
    % cluster AoD and AoA in az-el
    aod_az = aod_az_list(idx_clust);
    aod_el = aod_el_list(idx_clust);
    aoa_az = aoa_az_list(idx_clust);
    aoa_el = aoa_el_list(idx_clust);
    
    % azimuth spread
    tmp = (-spread_az:spread_az_res:spread_az);
    aod_list_az = tmp + aod_az;
    aoa_list_az = tmp + aoa_az;
    aod_list_az = aod_list_az(:) * pi/180;
    aoa_list_az = aoa_list_az(:) * pi/180;
    
    % elevation spread
    tmp = (-spread_el:spread_el_res:spread_el);
    aod_list_el = tmp + aod_el;
    aoa_list_el = tmp + aoa_el;
    aod_list_el = aod_list_el(:) * pi/180;
    aoa_list_el = aoa_list_el(:) * pi/180;
    
    % contributions of each ray
    for idx_ray_az_tx = 1:length(aod_list_az)
        for idx_ray_el_tx = 1:length(aod_list_el)
            for idx_ray_az_rx = 1:length(aoa_list_az)
                for idx_ray_el_rx = 1:length(aoa_list_el)
                    vtx = atx.get_array_response(aod_list_az(idx_ray_az_tx),aod_list_el(idx_ray_el_tx));
                    vrx = arx.get_array_response(aoa_list_az(idx_ray_az_rx),aoa_list_el(idx_ray_el_rx));
                    H = H + vrx * vtx';
                end
            end
        end
    end
end

% normalize channel energy
H = H ./ norm(H,'fro') .* sqrt(Nt*Nr);

% -------------------------------------------------------------------------
% 6. Compute mean.
% -------------------------------------------------------------------------
% coupling power of each beam pair over H (Gamma)
GG = abs(W' * H * F).^2;

% compute mean (mu)
mu = xi * 10 * log10(GG) + G_dB + EIRP_dBm - P_noise_dBm;

% -------------------------------------------------------------------------
% 7. Draw variance.
% -------------------------------------------------------------------------
% linear estimator
var = alpha * mu + beta;

% add random Gaussian noise
var = var + sqrt(nu_squared) .* randn(num_rx,num_tx);

% ensure variance is non-negative
var(var < 0) = 0;

% -------------------------------------------------------------------------
% 8. Log-normal realization of self-interference.
% -------------------------------------------------------------------------
INR = mu + sqrt(var) .* randn(num_rx,num_tx);

% -------------------------------------------------------------------------
% 9. Bound INR if desired.
% -------------------------------------------------------------------------
INR_bounds = [-Inf,Inf];
INR(INR < INR_bounds(1)) = INR_bounds(1);
INR(INR > INR_bounds(2)) = INR_bounds(2);

% -------------------------------------------------------------------------
% A. Plot resulting CDF of INR across beam pairs in the codebooks.
% -------------------------------------------------------------------------
[f,x] = ecdf(INR(:));
figure(1);
plot(x,f,'k-');
grid on;
grid minor;
xlabel('Self-Interference, INR (dB)');
ylabel('Cumulative Probability');
axis tight;

% -------------------------------------------------------------------------
% B. Plot resulting heatmap of INR across beam pairs in the codebooks.
% -------------------------------------------------------------------------
figure(2);
imagesc(INR);
xlabel('Transmit Beam Index');
ylabel('Receive Beam Index');
c = colorbar('EastOutside');
c.Label.Interpreter = 'latex';
c.Label.String = ['Self-Interference, INR (dB)'];
axis equal tight

%% -------------------------------------------------------------------------
% C. Plot only azimuth cut.
% -------------------------------------------------------------------------
idx_tx_az = find(tx_dir_az_el_deg(:,2) == 0);
idx_rx_az = find(rx_dir_az_el_deg(:,2) == 0);

figure(3);
surf(tx_dir_az,rx_dir_az,INR(idx_rx_az,idx_tx_az));
xlabel('Transmit Azimuth (deg.)');
ylabel('Receive Azimuth (deg.)');
c = colorbar('EastOutside');
c.Label.Interpreter = 'latex';
c.Label.String = ['Self-Interference, INR (dB)'];
axis equal tight
% shading interp;
view(0,90);