function [Avgdiff]=cdr_dereverc(i,u)

addpath(genpath('lib'));

%% filterbank initialization
cfg.K = 512; % FFT size
cfg.N = 10; % frame shift
cfg.Lp = 25; % prototype filter length
%p=IterLSDesign(cfg.Lp,cfg.K,cfg.N);
load('lib/filterbank/prototype_K512_N128_Lp1024.mat');

%% algorithm and scenario configuration
cfg.fs = 16000;      % sampling rate [Hz]
cfg.c = 342;         % speed of sound [m/s]
cfg.d_mic =  0.2;   % mic spacing [m]

% all estimators except estimate_cdr_nodoa require the TDOA of the signal; make sure
% to adapt this when loading another wave file

cfg.nr.lambda = 0.95; % smoothing factor for PSD estimation
cfg.nr.mu = 1.3;     % noise overestimation factor
cfg.nr.floor = 0.1;  % minimum gain
%cfg.nr.alpha = 1; cfg.nr.beta = 1; % power subtraction
cfg.nr.alpha = 2; cfg.nr.beta = 0.5; % magnitude subtraction
%cfg.nr.alpha = 2; cfg.nr.beta = 1; % Wiener filter

%cfg.estimator = @estimate_cdr_unbiased;           % unbiased estimator (CDRprop1)
%cfg.estimator = @estimate_cdr_robust_unbiased;    % unbiased, "robust" estimator (CDRprop2)
cfg.estimator = @estimate_cdr_nodoa;              % DOA-independent estimator (CDRprop3)
%cfg.estimator = @estimate_cdr_nodiffuse;          % noise coherence-independent estimator (CDRprop4; does not work for TDOA -> 0!)

%% preparation
% [i,fs_in] = audioread(filename_1);
% [u,fs_in] = audioread (filename_2);
x=[i,u];

%% Signal processing
% The algorithm itself is real-time capable, i.e., no processing of the entire
% utterance is necessary. Here however, for efficiency of the MATLAB implementation,
% the entire signal is processed at once.

fprintf('Performing signal enhancement... ');tic;

% analysis filterbank
X=DFTAnaRealEntireSignal(x,cfg.K,cfg.N,p);

% estimate PSD and coherence
Pxx = estimate_psd(X,cfg.nr.lambda);
Cxx = estimate_cpsd(X(:,:,1),X(:,:,2),cfg.nr.lambda)./sqrt(Pxx(:,:,1).*Pxx(:,:,2));

% frequency = (125:13:3453)'; % frequency axis
frequency = linspace(0,cfg.fs/2,cfg.K/2+1)'; % frequency axis
minf=min(frequency);
maxf=max(frequency);
% define coherence models
Css = exp(1j * 2 * pi * frequency );              % target signal coherence; not required for estimate_cdr_nodoa
Cnn = sinc(2 * frequency * cfg.d_mic/cfg.c); % diffuse noise coherence; not required for estimate_cdr_nodiffuse

% apply CDR estimator 
SNR = cfg.estimator(Cxx, Cnn, Css);
SNR = max(real(SNR),0);
Diff = 1./(SNR+1);
Avgdiff= (sum(sum(Diff)))/(257*(maxf-minf))

end