# Multichannel Acoustic Source Localization

This repository contains code to implement and evaluate various algorithms for localizing multiple acoustic sources using recordings from a microphone array. 
Please refer the [Report](https://github.com/gambiTarun/Multichannel-Acoustic-Source-Localization/blob/main/Report.pdf) for detailed information.

## Problem Statement

Given microphone recordings from an array, localize multiple acoustic sources by solving the inverse problem. 

## Methods Implemented

- Spectrum based beamforming
  - MVDR (Capon's) beamformer
  - Conventional (Bartlett) beamformer
- Deconvolution approaches
  - DAMAS
  - SC-DAMAS
  - Covariance matrix fitting
- Delay and sum beamforming

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

- Combination of methods needed for robust performance
- Extend to non-free field propagation environments  
- Incorporate 3D localization

## References

Refer the attached Report for references
