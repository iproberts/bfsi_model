# About

This repository contains code related to the following paper on modeling self-interference in full-duplex millimeter wave (mmWave) wireless communication systems.

[1] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Spatial and Statistical Modeling of Multi-Panel Millimeter Wave Self-Interference," Submitted to _IEEE J. Sel. Areas Commun._, 2023.

The work of [1] is based on our prior work in the following paper.

[2] I. P. Roberts, A. Chopra, T. Novlan, S. Vishwanath, and J. G. Andrews, "Beamformed Self-Interference Measurements at 28 GHz: Spatial Insights and Angular Spread," _IEEE Trans. Wireless Commun._, Nov. 2022.

Using the code in this repo, which is based on the model presented in [1], users can draw statistical realizations of self-interference in mmWave full-duplex systems. This can enable them to:
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

# What's The Difference Between [1] and [2]?

[1] and [2] can both be used to draw statistical realizations of self-interference in full-duplex mmWave systems. 

[1] is a more accurate in doing so in the sense that it captures spatial characteristics of self-interference whereas [2] is purely statistical. In other words, for a particular transmit and receive beam, [1] can produce realizations of self-interference that more closely align with what one would see in practice (based on our measurements). 

# Model Summary

Our statistical model of self-interference is based on two characteristics observed in our measurements:
1. On a large scale (at a high level), there is a connection between the steering directions of the transmit and receive beams and the degree of self-interference incurred. Broadly speaking, some transmit and receive directions tend to incur high self-interference while others tend to incur low self-interference.
2. On a small scale (within small spatial neighborhoods), the system incurs seemingly random amounts of self-interference. Slightly shifting the transmit and receive steering directions can dramatically alter the degree of self-interference coupled.

![A block diagram of our model.](https://user-images.githubusercontent.com/52005199/221431446-ae3a8393-2c4b-41e8-a66c-fcdf258f63e4.svg)

We leverage these large-scale and small-scale characteristics to construct a stochastic model of self-interference that both statistically and spatially aligns with our measurements. 
A block diagram summarizing our model is shown above.
For particular transmit and receive beams, a mean parameter is estimated, which dictates the location of the distribution from which self-interference is drawn.
The variance of this distribution is dictated by the mean parameter and other model parameters.
This approach allows our model to capture the large-scale spatial trends in self-interference along with the small-scale variability observed over small spatial neighborhoods.
With appropriate parameterization, our model has the potential to be extended to other systems and environments beyond our own. 

![Coupling clusters comprising the self-interference channel.](https://user-images.githubusercontent.com/52005199/221431451-9f7bec04-4659-4b9c-95d5-97deeb3f2345.svg)

To construct our model, we uncovered a coarse geometric approximation of the self-interference channel from within our measurements, which suggests that the dominant coupling between the transmit and receive arrays manifests as clusters of rays in a far-field manner (as illustrated above), rather than in a idealized near-field, spherical-wave fashion.
This is a novel finding that can steer future work aiming to model self-interference MIMO channels in full-duplex mmWave systems.

# Contents

This repo contains the following MATLAB code:
 - a main script illustrating example usage
 - an `array` object that can be used to create and interface with arbitrary antenna arrays

# Example Usage

The measurements and statistical characterization in [1] are particularly useful for drawing realizations of self-interference levels that a full-duplex mmWave platform may incur when transmitting while receiving in-band. We use interference-to-noise ratio (INR) to describe the level of self-interference experienced by such a system. When INR < 0 dB, the system is noise-limited, and when INR > 0 dB, it is self-interference-limited.

To draw a realization of mmWave self-interference (i.e., of INR) for a random transmit beam and random receive beam, one can execute the following in MATLAB using our code.

```
[m,s] = get_normal_params_min(0,0) % mean and variance
INR_dB = normrnd(m,sqrt(s)) % INR in dB
```

The first line fetches a mean `m` and variance `s`, which have been fitted from the measured INR distribution (see Fig. 5 in [1]).

The second line draws a realization of INR (in dB) from a normal distribution with mean `m` and standard deviation `sqrt(s)`.


# .mat Files

The `mat/` folder contains five `.mat` files corresponding to the five tables in our paper, each of which contains a lookup table (matrix) of fitted parameters. 

For instance, `params_normal_min.mat` corresponds to Table II in our paper and contains the following variables:
- `lut_mu_min` - a matrix where the `(i,j)`-th entry is the fitted mean of the `j`-th azimuthal neighborhood and `i`-th elevational neighborhood
- `lut_var_min` - a matrix where the `(i,j)`-th entry is the fitted variance of the `j`-th azimuthal neighborhood and `i`-th elevational neighborhood
- `delta_az_list` - the azimuthal neighborhood size of each row in `lut_mu_min` and `lut_var_min`
- `delta_el_list` - the elevational neighborhood size of each row in `lut_mu_min` and `lut_var_min`

# Questions and Feedback

Feel free to reach out to the corresponding author of [1] with any questions or feedback.

# Acknowledgments

This work has been supported by the National Science Foundation Graduate Research Fellowship Program (Grant No. DGE-1610403). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
