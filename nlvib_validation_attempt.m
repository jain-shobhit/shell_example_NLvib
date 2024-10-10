clearvars;
close all;
clc;
% This file attempts to perform frequeny response analysis via NLvib 1.4  on a 1320 degrees-of-freedom finite element model of a shallow-curved arch given in Example 6.3 of [1]
% [1] S. Jain & G. Haller, How to compute invariant manifolds and their reduced dynamics in high-dimensional finite-element models. Nonlinear Dyn. 107 (2022) 1417–1450. 
% Please install NLvib 1.4

addpath('../../SRC');
addpath('../../SRC/MechanicalSystems');
%% Build model
load('vonKarmanShell.mat','M','C','K','p','E','f_ext');


% Define oscillator as system with polynomial stiffness nonlinearities
system = ...
    System_with_PolynomialStiffnessNonlinearity(M,C,K,p,E,f_ext);

% Number of degrees of freedom
n = system.n;

% Frequency of first linear mode
[phi1,om12] = eigs(K,M,1,'sm'); % = 3.05 (plausible, cf. Fig. 6)
om1 = sqrt(om12);

%% Compute frequency response using harmonic balance

% Analysis parameters
analysis = 'FRF';
H = 10;      % harmonic order
N = 4*H+1;  % number of samples sufficient for max. cubic-degree polynomials, cf. Appendix A in https://doi.org/10.1016/j.ymssp.2019.106503
Om_s = 0.9*om1; % start frequency
Om_e = 1.1*om1; % end frequency

% Initial guess (solution of underlying linear system)
Q1 = (-Om_s^2*M + 1i*Om_s*C + K)\system.Fex1;
y0 = zeros((2*H+1)*length(Q1),1);
y0(length(Q1)+(1:2*length(Q1))) = [real(Q1);-imag(Q1)];
qscl = max(abs(y0));

% Solve and continue w.r.t. Om
ds = .05;
Sopt = struct('Dscale',[qscl*ones(size(y0));Om_s],'dynamicDscale',1);
tmp = tic;
X_HB = solve_and_continue(y0,...
    @(X) HB_residual(X,system,H,N,analysis),...
    Om_s,Om_e,ds,Sopt);
tFRC = toc(tmp);

