# MESD toolbox

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations.  By downloading and/or installing this software and associated files on your computing system you agree to use the software under the terms and condition as specified in the License agreement.

If this toolbox has been useful for you, please cite [1,2,3].

## Using the MESD toolbox

### About

This MATLAB toolbox implements the minimal expected switch duration (MESD) metric to evaluate the performance of neural decoders for auditory attention detection in the context of neuro-steered hearing prostheses as published in [1,2]. The ESD metric is an interpretable, single-number metric that combines both accuracy and decision window length. It allows easy comparison between neural decoders based on a interpretable criterion.

Developed and tested in MATLAB R2018b.

### Documentation

All functions are documented properly in their respective m-files. Additional documentation  and examples can be found in the [doc](doc/) folder, which contains a [manual](doc/manual.pdf) in pdf format and a [MESD demo file](doc/mesdDemo.m) to illustrate the usage of the various functions. A quick start guide is provided in the next section.
 
### Quick start guide
 
All functions needed to compute the MESD metric can be found in the [MESD folder](mesd-toolbox/). Before starting, make sure that this folder is added to the MATLAB path.

**Computation of the MESD performance metric**

The MESD is computed in four steps:

 1. Construction of the p(tau)-performance curve by interpolating through the evaluated (on real EEG and audio data) (tau_i,p_i)-points (decision window length accuracy).
 2. Optimization of the Markov chain in the number of states N for each sampled tau on the p(tau)-performance curve.
 3. Computation of the expected switch duration ESD(p(tau),tau,N_tau) per sampled tau and corresponding optimal number of states N_tau.
 4. The MESD is equal to the minimal expected switch duration over all evaluated ESD's.
 
These steps are implemented in the *main*-function [computeMESD.m](mesd-toolbox/computeMESD.m). Given the evaluated (tau_i,p_i)-performance points `(tau,p)`, the MESD can be computed with:

     mesd = computeMESD(tau,p);
 The default hyperparameter values are P_0 = 0.8 (confidence level), c = 0.65 (lower bound confidence interval), N_min = 5 (minimal number of states) and K = 1000 (number of samples evaluated on the performance curve). These hyperparameters can be adapted via extra arguments in the `computeMESD`-function.
 
 The toolbox also provides a [computeESD.m](mesd-toolbox/computeESD.m)-function to compute the ESD for a given performance pair (tau,p):
 
     esd = computeESD(tau,p);

**Designing an optimal Markov chain model for a neuro-steered hearing prosthesis** 

In Section *II-C*, a methodology is proposed to design an optimal Markov chain model for an adaptive gain control system in a neuro-steered hearing prosthesis. For a fixed accuracy p and hyperparameters P_0, c and N_min, the optimal number of states can be found with:

     Nopt = optimizeMarkovChain(p,Nmin,P0,c);


The optimal model for a certain neural decoder (represented by evaluated (tau_i,p_i)-points) can be identified by extra outputs of the `computeMESD`-function:

     [mesd,Nopt,tauOpt,pOpt] = computeMESD(tau,p,'Nmin',Nmin,'P0',P0,'c',c);

 ## References
 
[1] S. Geirnaert, T. Francart, and A. Bertrand, "An Interpretable Performance Metric for Auditory Attention Decoding Algorithms in a Context of Neuro-Steered Gain Control," IEEE Transactions on Neural Systems and Rehabilitation Engineering, 28(1), 307-317 (2020)

[2] S. Geirnaert, T. Francart, and A. Bertrand, "A New Metric to Evaluate Auditory Attention Detection Performance Based on a Markov Chain," Proc. European Signal Processing Conference (EUSIPCO), A Coruña, Spain, Sept. 2019

[3] S. Geirnaert, T. Francart, and A. Bertrand, "MESD toolbox," August 2019, Available online, URL: https://github.com/exporl/mesd-toolbox

