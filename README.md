# Multichannel Acoustic Source Localization

This project delves into the study of various methods to map acoustic sources using recordings from a microphone array. The longstanding problem of acoustic source localization is tackled using methodologies ranging from time-delay, beamforming, spectrum estimation, deconvolution, deep learning, and more. Acoustic sources are simulated using the solution for the Helmholtz wave equation. The project evaluates the robustness of implemented methods against variables like noise level, sensor configurations, wave propagation frequency, and more. Potential applications span from multi-source detection in HCI to locating neuron firing sources in the brain. 

Please refer the attached [Report](https://github.com/gambiTarun/Multichannel-Acoustic-Source-Localization/blob/main/Report.pdf) for detailed information.

## Problem Statement

Sound source localization originates from the human ability to pinpoint sound sources in 3D space using auditory data. The thesis implements and tests robustness of various DOA estimation methods mentioned in literature. Source locations are treated as an inverse problem using multichannel data, modeled on wave propagation rooted in Helmholtz equations.

## Methods Implemented

### Spectrum based beamforming

#### MVDR (Capon's) beamformer
The Minimum Variance Distortionless Response (MVDR) or Capon's beamformer optimizes the array response in a way that it has a unit gain at the direction of interest, while minimizing the output power. This results in higher resolution and better interference suppression compared to traditional beamformers.

#### Conventional (Bartlett) beamformer
The Conventional or Bartlett beamformer is a straightforward method where the outputs of each sensor are coherently added. This leads to a maximum response from the direction of the source of interest. It's widely used because of its simplicity, though it can be limited in resolution compared to other advanced beamforming techniques.

### Deconvolution approaches

#### DAMAS
The Deconvolution Approach for the Mapping of Acoustic Sources (DAMAS) algorithm is designed to iteratively deconvolve the Point Spread Function (PSF) from the beamformer output spectrum. This helps in extracting the actual source distribution from the delay-and-sum (DAS) beamformer spectrum, providing improved source localization especially in scenarios with multiple closely-spaced sources.

#### SC-DAMAS
The Sparse Component DAMAS (SC-DAMAS) is an enhancement to the traditional DAMAS, adding a regularization term that exploits the sparsity of the acoustic sources in the spatial domain. This addition helps in achieving faster convergence and potentially better source resolution.

#### Covariance matrix fitting
This method involves formulating an optimization problem to fit the measured covariance matrix with a model covariance matrix. This matrix is constructed from potential source locations, enabling the localization of acoustic sources even in environments with high levels of noise.

### Delay and sum beamforming
Delay-and-sum is one of the most basic forms of beamforming. Signals from each sensor are delayed according to the estimated time-of-arrival of the wavefront from the direction of the source, and then summed together. The method is straightforward and computationally effective, but can face challenges in environments with multiple sources or noise.

## Experimental Setups

- 2D propagation with circular microphone array
- 3D propagation with planar microphone array  

## Evaluation

Methods evaluated on robustness w.r.t:

- Additive noise levels
- Number of microphones 
- Acoustic frequency
- Microphone array configuration
- Distance between sources and microphones

## Results

- MVDR beamforming effective for well-separated sources
- Deconvolution approaches better for close sources
- Covariance matrix fitting best at high noise levels
- No single technique uniformly better across scenarios

## Conclusion

Our study encompassed various algorithms for source localization, ranging from spectrum-based beamforming approaches to deconvolution methods, and direct optimization using Covariance Matrix Fitting (CMF). Key observations include:

- MVDR excels particularly in low noise environments, with more sensors or at higher frequencies. Its performance is influenced by the sensor array shape, especially in Setup 1 2.1.
- While MVDR offers better resolution when sources are well separated, Bartlett beamforming presents as a simpler alternative with marginally lower resolution.
- Deconvolution methods like DAMAS and SC-DAMAS, built atop the delay-and-sum beamformer, aim to enhance resolution by deconvolving the spectrum. Notably, SC-DAMAS outperforms other methods when sources are further from sensors.
- In noisier environments, CMF surpasses DAMAS and SC-DAMAS in performance.
- Future research can delve into non-free-field scenarios, more representative of real-world conditions, though this introduces challenges in estimating medium-specific acoustic properties.

## References

Refer the attached Report for references
