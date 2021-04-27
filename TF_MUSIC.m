clc
clear
clear all

tic
fs      = 1;                    % Sampling frequency
Samples       = 512;            % Number of samples
FreqBins       = 512;           % Number of frequency bins
perc    = 0.8   ;               % Percentage of the STFD maximum power to select high-energy (t,f) points
theta   = -90:0.1:90;    
doa = [-30 -25 40 80]/180*pi;   %Direction of arrival in radians
M = 6;                          %Number of array elements
P = 4;                          %The number of signal
lambda = 100;                   %Wavelength
d = lambda/2;                   %Element spacing
SNR = 20    ;                   %Signal-noise-ratio
A = zeros(P,M);                 %To creat a matrix with P row and M column
for k = 1:P
    A(k,:) = exp(-j*2*pi*d*sin(doa(k))*(0:M-1)/lambda); %Steering matrix
end
A=A';
f_init  = [0.0 0.5 0.0 0.4];       % LFM initial frequencies of P sources
f_end   = [ 0.5 0.0 0.5 0.1];
SS = signal_model(P, f_end, f_init, Samples, fs); % n LFM signal generation
S = A*SS;
S = S + awgn(S,SNR);            %Insert Gaussian white noise
Rx = S*S';                      %Data covariance matrix
Rxx = S*S'/(Samples);           %Data covariance matrix

%% Multisensor Time-Frequency Distribution (MTFD)
        D = cell(M,M);
        for i = 1:M
            for j = 1:M
                s1 = S(i,:);        % signal from one microphone
                s2 = S(j,:);        % signal from the next microphone
                
                N = length(s1);     
                if(~mod(N,2))
                    L = N-1;        % Length of the lag window (should be odd)
                else
                    L = N;
                end
                MM = 2^nextpow2(N);  % new length of a signal divisible by 2
                PP = nextpow2(N);    % 
                K = 2*(2^PP);        % doubling the length of a signal for FFT calculations
                ww1 = {'bartlett','barthannwin','blackman','blackmanharris','bohmanwin','chebwin','flattopwin','hamming',...
    'hann','nuttallwin','parzenwin','rectwin','taylorwin','triang','gausswin'}';
                win = 'blackmanharris';
                z1 = fft(s1,K);           % Transforming S1 into the frequency domain
                z1(2:K/2) = 2*z1(2:K/2);  % Accepting and doubling half of the frequencies
                z1(K/2+1:K) = 0;          % Rejecting half of the frequencies
                z1 = ifft(z1);            % Transforming Z into the time domain
                z2 = fft(s2,K);           % Transforming S2 into the frequency domain
                z2(2:K/2) = 2*z2(2:K/2);  % Accepting and doubling half of the frequencies
                z2(K/2+1:K) = 0;          % Rejecting half of the frequencies
                z2 = ifft(z2);            % Transforming Z into the time domain

                K_TL = zeros(K, N);
                L_half = fix(L/2);                
                if(sum(strcmp(win,ww1)))
                    g = window(str2func(win),L);
                elseif(sum(strcmp(win,ww2)))
                    g = window(str2func(win),L,OPT);
                end
                
                for n = 1:K
                    for tau = -L_half:L_half
                        G  = g(1 + tau + L_half);  
                        Z1 = z1(1 + rem(2*K + n - 1 + tau, K));
                        Z2 = z2(1 + rem(2*K + n - 1 - tau, K));
                        mm = 1 + rem(2  *K + tau, MM);
                        K_TL(mm,n) = G*Z1*conj(Z2);
                    end
                end
                TFR  = zeros(MM, N);
                temp = fft(K_TL,MM);
                TFR(:,2:N) = temp(:,1:N-1);
                D{i,j} = TFR(:,2:N);
            end      
        end
       
        %% Averaged Auto-TFD
        D_avg = zeros(FreqBins,Samples-1);
        for ii = 1:M
              %Dmm = size(D{ii,ii})
              %Davg = size(D_avg)
            D_avg = D{ii,ii}+D_avg; 
        end
        D_avg = D_avg./M;
        
        %% Selection of high-energy (t,f) points
        thr = perc*max(max(D_avg));
        Tr = abs(D_avg) >= thr;
        [F_trace, ~] = find(Tr);
        n_p = length(F_trace);
        D_s = zeros(M, M);
        for m1 = 1:M
            for m2 = 1:M
                D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
            end
        end
%% TF-MUSIC Algorithm

theta   = theta/180*pi;
theta_N = length(theta);
[ee,~]  = eig(D_s);        % Find the eigenvalues and eigenvectors of STFD
UN_TF   = ee(:,P+1:M);    % Estimate/Selection of noise subspace

P = zeros(1,theta_N);
for ii = 1:theta_N
    a_theta = exp(-1j*2*pi*(d/lambda)*sin(theta(ii))*(0:M-1));
    PP_TF   = (a_theta*UN_TF)*(UN_TF'*a_theta');
    P(ii)   = abs(1/PP_TF);
end
PTFMUSIC = 10*log10( P/max(P));
figure
plot(theta/pi*180,PTFMUSIC,'LineWidth',1.5)
xlabel('DoA angle, \theta ','FontSize',14)
ylabel('P, dB','FontSize',14)
title('TF-MUSIC Algorithm','FontSize',16)
grid on
grid minor

%%
% for i = 1:length(DOA_MTFD)
%     e1 = -30 - DOA_MTFD(1,i);
%     e2 = -25 - DOA_MTFD(1,i);
%     e3 = 40  - DOA_MTFD(1,i);
%     e4 = 80  - DOA_MTFD(1,i);

%     minE = min(e1,e2,e3,e4)
%     if(-30)
%     else if
% end
% rmse=sqrt(mean((a(:)-b(:)).^2));
