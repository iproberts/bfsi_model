# About

This repository contains code related to the following paper on modeling self-interference in full-duplex millimeter wave (mmWave) wireless communication systems.

[1] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Spatial and Statistical Modeling of Multi-Panel Millimeter Wave Self-Interference," Submitted to _IEEE J. Sel. Areas Commun._, 2023, [PDF](https://ianproberts.com/pdf/pub/simodel.pdf).

The work of [1] is based on measurements collected in our prior work in the following paper.

[2] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Beamformed Self-Interference Measurements at 28 GHz: Spatial Insights and Angular Spread," _IEEE Trans. Wireless Commun._, Nov. 2022, [PDF](https://ianproberts.com/pdf/pub/bfsi.pdf).

Using the code in this repo, which is based on the model presented in [1], users can draw statistical realizations of self-interference in mmWave full-duplex systems. This can allow them to:
 - conduct statistical analyses of full-duplex mmWave communication systems;
 - develop methods to mitigate self-interference in full-duplex mmWave communication systems;
 - evaluate solutions for full-duplex mmWave communication systems.

This work is based on nearly 6.5 million measurements of self-interference taken at 28 GHz in an anechoic chamber using two colocated 256-element phased arrays mounted on separate sides of an equilateral triangular platform. Please see [1] and [2] for details for a summary of these measurements.

If you use this code or our paper in your work, please cite [1] with the following BibTeX.

```
@ARTICLE{roberts_si_model_2023,
    author={I. P. Roberts and A. Chopra and T. Novlan and S. Vishwanath and J. G. Andrews},
    journal={Submitted to IEEE J.~Sel.~Areas~Commun.},
    title={Spatial and Statistical Modeling of Multi-Panel Millimeter Wave Self-Interference}, 
    year=2023,
}
```

Related work can be found at https://ianproberts.com.

# What is Self-Interference? 

When a transceiver (a wireless device) attempts to transmit and receive at the same time using the same frequency spectrum, some of its transmitted signal will leak into its receiver, corrupting reception of a desired signal.
This undesired leakage (or coupling) is called "self-interference". 
Transmit and receiving at the same time using the same frequency spectrum is called "full-duplex" operation.

This work is particularly focused on self-interference in full-duplex systems operating at mmWave frequencies (roughly 30 GHz to 100 GHz).
In mmWave systems, dense antenna arrays containing dozens or even hundreds of individual antennas are used to overcome high path loss at these high carrier frequencies. 
They do this by forming very directional beams, focusing energy in a particular direction to increase received signal power. 

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/222208767-c359fd9d-0fe0-4814-a56b-46d1f3fd306d.svg"/>
</p>

This work is therefore interested in how much self-interference is coupled in full-duplex mmWave systems when using particular transmit and receive beams. Some transmit and receive beams will couple higher self-interference than others; this depends on the steering direction of the beams and on the (unknown) underlying self-interference channel. This work aims to model self-interference statistically and spatially in order to draw realizations that align with the levels that one would see in practice, based on our measurements.

# What's The Difference Between [1] and [2]?

[1] and [2] can both be used to draw statistical realizations of self-interference in full-duplex mmWave systems. 

[1] is a more accurate in doing so in the sense that it captures spatial characteristics of self-interference whereas [2] is purely statistical. In other words, for a particular transmit and receive beam, [1] can produce realizations of self-interference that more closely align with what one would see in practice (based on our measurements). 

# Model Summary

Our statistical model of self-interference is based on two characteristics observed in our measurements:
1. On a large scale (at a high level), there is a connection between the steering directions of the transmit and receive beams and the degree of self-interference incurred. Broadly speaking, some transmit and receive directions tend to incur high self-interference while others tend to incur low self-interference.
2. On a small scale (within small spatial neighborhoods), the system incurs seemingly random amounts of self-interference. Slightly shifting the transmit and receive steering directions can dramatically alter the degree of self-interference coupled.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/221431446-ae3a8393-2c4b-41e8-a66c-fcdf258f63e4.svg"/>
</p>

We leverage these large-scale and small-scale characteristics to construct a stochastic model of self-interference that both statistically and spatially aligns with our measurements. 
A block diagram summarizing our model is shown above.
For particular transmit and receive beams, a mean parameter is estimated, which dictates the location of the distribution from which self-interference is drawn.
The variance of this distribution is dictated by the mean parameter and other model parameters.
This approach allows our model to capture the large-scale spatial trends in self-interference along with the small-scale variability observed over small spatial neighborhoods.
With appropriate parameterization, our model has the potential to be extended to other systems and environments beyond our own. 

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/221431451-9f7bec04-4659-4b9c-95d5-97deeb3f2345.svg"/>
</p>

To construct our model, we uncovered a coarse geometric approximation of the self-interference channel from within our measurements, which suggests that the dominant coupling between the transmit and receive arrays manifests as clusters of rays in a far-field manner (as illustrated above), rather than in a idealized near-field, spherical-wave fashion.
This is a novel finding that can steer future work aiming to model self-interference MIMO channels in full-duplex mmWave systems.

# Contents

This repo contains the following MATLAB code:
 - a main script `main.m` illustrating example usage
 - an `array` object that can be used to create and interface with arbitrary antenna arrays

# Example Usage

Suppose a full-duplex mmWave base station employs a codebook of transmit beams that it uses to serve downlink and a codebook of receive beams that it uses to serve uplink. 

The degree of self-interference incurred by the base station will depend on the transmit beam and receive beam that it uses to serve downlink and uplink users. Each transmit-receive beam pair will couple a unique amount of self-interference. 

We can draw a statistical realization of this coupling for each beam pair across the codebooks using the script `main.m`, which implements our model in MATLAB.

### Set System and Model Parameters

The first step is to set the desired system and model parameters. By default, the parameters provided in [1, Table II] can be used. 

```
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
```

### Set Transmit and Receive Arrays

The transmit and receive array sizes can be modified by changing the following lines. Our arrays were 16x16 planar arrays, for instance. It is not clear how well our model will generalize to other array sizes.

```
% planar array size
M = 16; % number of rows
N = 16; % number of columns
```

### Define Transmit and Receive Codebooks

The transmit and receive codebooks used at the full-duplex base station can be set by defining the steering directions of the codebooks' beams. The steering direction of each beam contains a component in azimuth and elevation (see [1, Fig. 2]). 

In the example below, the transmit and receive codebooks are identical. Each codebook has beams that span in azimuth from -56 deg. to 56 deg. in 8 deg. steps and in elevation from -8 deg. to 8 deg. in 8 deg. steps. This amounts to a total of 45 beams in each codebook, meaning there are 2025 transmit-receive beam pairs.

```
% transmit directions (e.g., transmit codebook)
tx_dir_az = flip([-56:8:56]); % flips are not necessary, just for plotting convenience
tx_dir_el = flip([-8:8:8]);

% receive directions (e.g., receive codebook)
rx_dir_az = tx_dir_az; % assume same for simplicity
rx_dir_el = tx_dir_el;
```

The azimuth-elevation of each transmit and receive steering direction can then be populated as follows.

```
% each transmit az-el pair
tx_dir_az_el_deg = [repelem(tx_dir_az(:),length(tx_dir_el)),repmat(tx_dir_el(:),length(tx_dir_az),1)];
num_tx = length(tx_dir_az_el_deg(:,1));

% each receive az-el pair
rx_dir_az_el_deg = [repelem(rx_dir_az(:),length(rx_dir_el)),repmat(rx_dir_el(:),length(rx_dir_az),1)];
num_rx = length(rx_dir_az_el_deg(:,1));
```

To steer in a particular direction, we use conjugate beamforming (or matched filter beamforming). This can be achieved via the following, which uses the provided `array` object to fetch the array response vector of each array in each steering direction.

```
% array response vectors in each steering direction
Atx = atx.get_array_response(tx_dir_az_el_deg(:,1)*pi/180,tx_dir_az_el_deg(:,2)*pi/180);
Arx = arx.get_array_response(rx_dir_az_el_deg(:,1)*pi/180,rx_dir_az_el_deg(:,2)*pi/180);

% transmit and receive codebooks use conjugate beamforming
F = Atx;
W = Arx;
```

Here, the `i`-th column in the matrix `F` contains the transmit beamforming vector that steers in the `i`-th direction of `tx_dir_az_el_deg`.

It is very important that the energy of the transmit and receive beams be normalized properly. To ensure this is the case, the following is used.

```
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
```

### Constructing the Clustered Self-Interference Channel Approximation

The following portion of the code constructs the channel matrix used to estimate the neighborhood mean parameter.

```
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
```

### Computing the Neighborhood Mean and Variance

The neighborhood mean for each transmit-receive beam pair can then be computed directly using the following.

```
% coupling power of each beam pair over H (Gamma)
GG = abs(W' * H * F).^2;

% compute mean (mu)
mu = xi * 10 * log10(GG) + G_dB + EIRP_dBm - P_noise_dBm;
```

Likewise, the neighborhood variance can be drawn as follows.

```
% linear estimator
var = alpha * mu + beta;

% add random Gaussian noise
var = var + sqrt(nu_squared) .* randn(num_rx,num_tx);

% ensure variance is non-negative
var(var < 0) = 0;
```

### Drawing Realizations of Self-Interference

Self-interference for each beam pair is then straightforwardly drawn as follows. Note that `INR` is the interference-to-noise ratio (INR) of self-interference and has units of decibels (dB).

```
INR = mu + sqrt(var) .* randn(num_rx,num_tx);
```

If desired, `INR` can be bounded by setting `INR_bounds`. By default it is not bounded.

```
INR_bounds = [-Inf,Inf];
INR(INR < INR_bounds(1)) = INR_bounds(1);
INR(INR > INR_bounds(2)) = INR_bounds(2);
```

### Results

The INR for each transmit-receive beam pair can then be plotted using the following.

```
figure(1);
imagesc(INR);
xlabel('Transmit Beam Index');
ylabel('Receive Beam Index');
c = colorbar('EastOutside');
c.Label.Interpreter = 'latex';
c.Label.String = ['Self-Interference, INR (dB)'];
axis equal tight
```

This will produce a figure similar to the following. Each pixel represents the degree of self-interference incurred by the corresponding transmit-receive beam pair.

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/221433216-4fd2667e-8ecc-4cc6-9eaa-01de70647421.svg"/>
</p>

The distribution of this can be plotted as follows, producing a plot similar to the one shown below.

```
[f,x] = ecdf(INR(:));
figure(2);
plot(x,f,'k-');
grid on;
grid minor;
xlabel('Self-Interference, INR (dB)');
ylabel('Cumulative Probability');
axis tight;
```

<p align="center">
<img src="https://user-images.githubusercontent.com/52005199/221433220-befc149f-68c3-4978-8fd0-c0edc9ff64ed.svg"/>
</p>

# Questions and Feedback

Feel free to reach out to the corresponding author of [1] with any questions or feedback.

# Acknowledgments

This work has been supported by the National Science Foundation Graduate Research Fellowship Program (Grant No. DGE-1610403). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
