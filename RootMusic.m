clear
clear all
clc

fs       = 1;                % Sampling frequency
Samples  = 512;              % Number of samples
FreqBins = 512;              % Number of frequency bins
perc     = 0.5;              % Percentage of the STFD maximum power to select high-energy (t,f) points
theta    = -90:0.5:90;
    
doa = [-30 -25 40 80]/180*pi;     %Direction of arrival in radians
w = [pi/2 pi/4 pi/3  pi/5]';      %Frequency (angular)

N = 100;                         %Snapshots
M = 6;                           %Number of array elements
P = 4;                           %The number of signal
lambda = 100;                    %Wavelength
d = lambda/2;                    %Element spacing
SNR = -15;                       %Signal-noise-ratio
A = zeros(P,M);                  %To creat a matrix with P row and M column
for k = 1:P
    A(k,:) = exp(-j*2*pi*d*sin(doa(k))*(0:M-1)/lambda); %Steering matrix
end
A=A';

Samples = 512;
fs=1;
f_init  = [0.0, 0.2,0.3,0.5];       % LFM initial frequencies of n sources
f_end   = [0.2, 0.0, 0.5,0.3]; 
SS = signal_model(P, f_end, f_init, Samples, fs); % n LFM signal generation
% SS = 10*sin(w*(1:N));
S = A*SS;

S = S + awgn(S,SNR);      
Rx = S*S'/(Samples);          %Data covariance matrix

Built_in_Root_MUSIC = rootmusicdoa(Rx,4)  % Built-in Root MUSIC algorithm results

[vv,~]  = eig(Rx);          % Find the eigenvalues and eigenvectors of STFD
NN   = vv(:,1:M-P);         % Estimate/Selection of noise subspac

En = NN*NN';    
b = zeros(2*(M - 1), 1);
for i = -(M-1):1:M-1
    b(i+M) = sum(diag(En, i));
end
b = flipud(b);
rts = roots(b);
distance = 1 - abs(rts);

for ii = 1: length(distance)
    for jj = 2 : ii
       if(abs(distance(jj-1))>=abs(distance(jj)))
         temp1 = distance(jj-1);
         temp2 = rts(jj-1);
         distance(jj-1) = distance(jj);
         distance(jj) = temp1; 
         rts(jj-1) = rts(jj);
         rts(jj) = temp2;        
       end
    end
end
DOA = zeros(1,P*2);
for n = 1 : P*2
     DOA(n) = asin((angle(rts(n))*lambda)/(2*pi*d))*180/pi;
end
Built_in_Root_MUSIC = DOA



