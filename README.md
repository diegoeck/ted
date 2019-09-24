# ted

The Experiment Design (TED)

## Description

This is a Matlab toolbox for *experiment design* for system identification of dynamic systems.
The toolbox solves implement two different parameterizations for the input signal:
* Discrete Spectra (realization as sum of sines)
* Continuous Spectra (realization as filtered white noise)

In order to ensure the obtained result is a Spectra ( Phi(w)>0 ), some constrains are imposed to the optimization problem. Both continuous with finite parameterization and Partial Correlation Approach are implemented.



## Install

Just add the folder ./src/ to the *path* of Matlab.

## Use

Please check the *examples* folder.

#### Example of Discrete Spectra

```matlab
% Define System
G = zpk([0.7],[0.8 0.9],1,1);
H = tf([1],[1],1);
sigma = 1;
M0 = tf2model(G,H,sigma);

% Define Energy, Length and Frequencies of Input Signal
energy = 1;
length = 1000;
Ini = 3;
Freq = 30;
W1 = [];
for a=1:Freq
    W1(a) = 10^(((a-1)*Ini/Freq)-Ini)*pi;
end

% Compute Spectra
[BkP] = experiment_Fu(M0,W1);
[R0] = experiment_R0_cos_ss(M0);
[r, invP, topt] = experiment_e(BkP, R0, energy)

% Realize Spectrum as Sum os Sines
u = experiment_realization(r, W1, length);
```


## Contributors

Diego Eckhard - diegoeck@ufrgs.br - @diegoeck
