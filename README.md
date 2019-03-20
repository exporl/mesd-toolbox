
# EST toolbox

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations.  By downloading and/or installing this software and associated files on your computing system you agree to use the software under the terms and condition as specified in the License agreement.

If this toolbox has been useful for you, please cite [1,2,3].

## Using the EST toolbox

### About

This MATLAB toolbox implements the Markov transit time (EST) metric to evaluate the performance of neural decoders for auditory attention detection in the context of neuro-steered hearing prostheses as published in [1,2]. The EST metric is an interpretable, single-number metric that combines both accuracy and decision time. It allows easy comparison between neural decoders based on a relevant criterion, independent of the evaluated window lengths.

Developed and tested in MATLAB R2018b.

### Documentation

All functions are documented properly in their respective m-files. Additional documentation  and examples can be found in the [doc](doc/) folder, which contains a  [manual](doc/manual.pdf) in pdf format and a [EST demo file](doc/estDemo.m) to illustrate  the usage of the various functions. A quick start guide is provided in the next section.
 
### Quick start guide
 
All functions needed to compute the EST metric can be found in the [EST folder](est-toolbox/). Before starting, make sure that this folder is added to the MATLAB path.

**Computation of the EST performance metric**

The EST is computed in four steps:

 1. Construction of the p(tau)-performance curve by interpolating through the evaluated (on real EEG and audio data) (tau_i,p_i)-points  (decision time, accuracy).
 2. Optimization of the Markov chain in the number of states N for each sampled tau on the p(tau)-performance curve.
 3. Computation of the transit time T(p(tau),tau,N_tau) per sampled tau and corresponding optimal number of states N_tau.
 4. The EST is equal to the minimal transit time over all evaluated transit times.
 
These steps are implemented in the *main*-function [computeEST.m](est-toolbox/computeEST.m). Given the evaluated (tau_i,p_i)-performance points `(tau,p)`, the EST can be computed with:

     est = computeEST(tau,p);
 The default hyperparameter values are P_0 = 0.9 (confidence level), c = 0.65 (lower bound confidence interval), N_min = 5 (minimal number of states) and K = 1000 (number of samples evaluated on the performance curve). These hyperparameters can be adapted via extra arguments in the `computeEST`-function.
 
**Designing an optimal Markov chain model for a neuro-steered hearing prosthesis** 

In Section *II.D*, a methodology is proposed to design an optimal Markov chain model for an adaptive gain control system in a neuro-steered hearing prosthesis. For a fixed accuracy p and hyperparameters P_0, c and N_min, the optimal number of states can be found with:

     Nopt = optimizeMarkovChain(p,Nmin,P0,c);


The optimal model for a certain neural decoder (represented by evaluated (tau_i,p_i)-points) can be identified by extra outputs of the `computeEST`-function:

     [est,Nopt,tauOpt,pOpt] = computeEST(tau,p,'Nmin',Nmin,'P0',P0,'c',c);

 ## References
 
[1] S. Geirnaert, T. Francart, and A. Bertrand,  "Expected Switching Time: an Interpretable Performance Metric to Evaluate Neural Decoders for Auditory Attention Detection," March 2019, *Internal Report*

[2] S. Geirnaert, T. Francart, and A. Bertrand,  "A New Metric to Evaluate Auditory Attention Detection Performance Based on a Markov Chain," March 2019, *Internal Report*

[3] S. Geirnaert, T. Francart, and A. Bertrand, “EST toolbox,” March 2019, Available online, URL: https://github.com/exporl/est-toolbox
