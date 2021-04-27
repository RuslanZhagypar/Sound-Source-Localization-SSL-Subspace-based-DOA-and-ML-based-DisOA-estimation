clc
clear all
tic
    
doa = [-30 -25 40 80]/180*pi;    %Direction of arrival in radians
N = 100;                         %Snapshots
w = [pi/2 pi/4 pi/3 pi/5]';      %Frequency (angular)
M = 6;                           %Number of array elements
P = length(doa);                 %The number of signal
lambda = 100;                    %Wavelength
d = lambda/2;                    %Element spacing
SNR = -15;                       %Signal-noise-ratio
A = zeros(P,M);                  %To creat a matrix with P row and M column
for k = 1:P
    A(k,:) = exp(-j*2*pi*d*sin(doa(k))*(0:M-1)/lambda); %Steering matrix
end
A=A';
%SS = 10*rand(1,1)*sin(100*rand(1,1)*w*(1:N));     %Simulated signal

Samples = 512;
fs=1;
f_init  = [0.0, 0.2,0.3,0.5];       % LFM initial frequencies of n sources
f_end   = [0.2, 0.0, 0.5,0.3];  

SS = signal_model(P, f_end, f_init, Samples, fs); % n LFM signal generation
% SS = 10*sin(w*(1:N));
S = A*SS;

% N = 0.1*rand(1,length(S));
% b=[1 0.5];
% a= [1 -0.2 +0.5 1.6];
% Nfiltered = filter(b,a,N);
% S = S + Nfiltered;         %Insert Gaussian colored noise

S = S + awgn(S,SNR);        %Insert Gaussian white noise
Rx = S*S';                  %Data covariance matrix

% J = flip(eye(M));         % Spatial Smoothing 
% Ry = J*conj(Rx)*J;
% Rx = Rx + Ry;
[N,V] = eig(Rx);            %Find the eigenvalues and eigenvectors of R

NN = N(:,1:M-length(doa));  %Estimate noise subspace  
theta = -90:0.5:90;         %Peak search
for ii = 1:length(theta)
    SS = zeros(1,length(M));
    for jj = 0 : M-1
        SS(1+jj) = exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
    end
    PP = (SS*NN)*(NN'*SS');
    PMUSIC(ii) = abs(1/PP);
end

PMUSIC = 10*log10(PMUSIC/max(PMUSIC)); %Spatial spectrum function
plot(theta,PMUSIC,'LineWidth',1.5)
xlim([-90 90])
%ylim([-50 1])
xlabel('DoA angle, \theta ','FontSize',14)
ylabel('P, dB','FontSize',14)
title('MUSIC Algorithm with White Noise','FontSize',16)
grid on
grid minor
toc