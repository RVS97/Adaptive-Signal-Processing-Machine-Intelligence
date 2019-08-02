%% EX 1.1


%% EX 1.2
% 1.2 a)
load sunspot.dat

x = sunspot(:, 2);
xMean = x - mean(x);
xDetrend = detrend(x);
logX = log(x+eps);
logXmean = logX - mean(logX);

[PSDmean,f1] = periodogram(xMean,hamming(length(xMean)),[],1);
[PSDdetrend,f2] = periodogram(xDetrend,hamming(length(xDetrend)),[],1);
[PSDlog,f3] = periodogram(logXmean,hamming(length(logXmean)),[],1);

figure
plot(f1,10*log10(PSDmean), 'linewidth', 2)
hold on
plot(f2,10*log10(PSDdetrend), 'linewidth', 2)
plot(f3,10*log10(PSDlog), 'linewidth', 2)
hold off
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('mean','detrend','log')
title('Periodogram Sunspot Data Series', 'Fontsize', 35)

% 1.2 b)
load('.\EEG_Data\EEG_Data_Assignment1.mat')

[PSDPOz,f4] = periodogram(POz,hamming(length(POz)),5*fs,fs);
figure
plot(f4,10*log10(PSDPOz), 'linewidth', 2)
xlabel('Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
xlim([0 60])
title('Periodogram for EEG data', 'Fontsize', 35)


POz1 = reshape(POz-mean(POz), fs, []);
POz5 = reshape(POz-mean(POz), fs*5, []);
POz10 = reshape(POz-mean(POz), fs*10, []);
PSDPOz1 = mean(periodogram(POz1,hamming(length(POz1)),5*fs,fs),2);
PSDPOz5 = mean(periodogram(POz5,hamming(length(POz5)),5*fs,fs),2);
PSDPOz10 = mean(periodogram(POz10,hamming(length(POz10)),5*fs,fs),2);

% TODO frequencies axis
figure
plot(f4,10*log10(PSDPOz1), 'linewidth', 2);
hold on
plot(f4,10*log10(PSDPOz5), 'linewidth', 2);
plot(f4,10*log10(PSDPOz10), 'linewidth', 2);
xlim([0 60]);
xlabel('Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('1s window','5s window','10s window','Location','northwest')
title('Periodogram for windowed EEG data', 'Fontsize', 35)

%% EX 1.3
% a)
wgn = randn(1,2048);
noiseSin = sin(2*pi*(1:2048)*0.4) + wgn;
noiseFilt = filter([1/4 1/4 1/4 1/4], 1, wgn); %change ?

[wgnB, lag1] = xcorr(wgn, 'biased');
wgnUnb = xcorr(wgn, 'unbiased');
[sinB, lag2] = xcorr(noiseSin, 'biased');
sinUnb = xcorr(noiseSin, 'unbiased');
[filtB, lag3] = xcorr(noiseFilt, 'biased');
filtUnb = xcorr(noiseFilt, 'unbiased');

PSDwB = fftshift(real(fft(ifftshift(wgnB))));
PSDwU = fftshift(real(fft(ifftshift(wgnUnb))));
PSDsB = fftshift(real(fft(ifftshift(sinB))));
PSDsU = fftshift(real(fft(ifftshift(sinUnb))));
PSDfB = fftshift(real(fft(ifftshift(filtB))));
PSDfU = fftshift(real(fft(ifftshift(filtUnb))));

% Use same lag for all?
figure
plot(lag1, wgnUnb, 'linewidth', 2);
hold on
plot(lag1, wgnB, 'linewidth', 2);
xlabel('Lag')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Unbiased','Biased')
title('ACF for WGN', 'Fontsize', 35)
figure
plot(lag2, sinUnb, 'linewidth', 2);
hold on
plot(lag2, sinB, 'linewidth', 2);
xlabel('Lag')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Unbiased','Biased')
title('ACF for noisy Sinewave', 'Fontsize', 35)
figure
plot(lag3, filtUnb, 'linewidth', 2);
hold on
plot(lag3, filtB, 'linewidth', 2);
xlabel('Lag')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Unbiased','Biased')
title('ACF for filtered Noise', 'Fontsize', 35)

normf = linspace(-1,1, length(lag1));
figure
plot(normf,PSDwU, 'linewidth', 2);
hold on
plot(normf,PSDwB, 'linewidth', 2);
xlabel('Normalised Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Unbiased','Biased')
title('PSD for WGN', 'Fontsize', 35)
figure
plot(normf,PSDsU, 'linewidth', 2);
hold on
plot(normf,PSDsB, 'linewidth', 2);
xlabel('Normalised Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Unbiased','Biased')
title('PSD for noisy Sinewave', 'Fontsize', 35)
figure
plot(normf,PSDfU, 'linewidth', 2);
hold on
plot(normf,PSDfB, 'linewidth', 2);
xlabel('Normalised Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Unbiased','Biased')
title('PSD for filtered noise', 'Fontsize', 35)

% b) take average and std of periodogram
N = 64;

wgn100 = randn(1000,N);
noiseSin100 = repmat(sin(2*pi*(1:N)*0.4) + sin(2*pi*(1:N)*0.15), 1000,1) + wgn100;

sinB100 = zeros(1000,2*N-1);

for i=1:1000
    sinB100(i,:) = xcorr(noiseSin100(i,:), 'biased');
end

PSDs100 = fftshift(real(fft(ifftshift(sinB100,2)')));

meanPSDs100 = mean(PSDs100,2);
stdPSDs100 = std(PSDs100,0,2);

normf = linspace(-0.5,0.5,length(meanPSDs100));
figure
plot(normf, PSDs100, 'linewidth', 2,'color','cyan');
hold on
plot(normf, meanPSDs100,'b', 'linewidth', 2);
xlabel('Normalised Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('PSD estimates and mean of noisy sinewave', 'Fontsize', 35)
figure
plot(normf, stdPSDs100, 'linewidth', 2);
xlabel('Normalised Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('STD of noisy sinewave PSD', 'Fontsize', 35)

% c : in dB

PSDs100dB = pow2db(PSDs100);
meanPSDs100dB = pow2db(meanPSDs100);
stdPSDs100dB = pow2db(stdPSDs100);

figure
plot(normf, PSDs100dB, 'linewidth', 2,'color','cyan');
hold on
plot(normf, meanPSDs100dB,'b', 'linewidth', 2);
xlabel('Normalised Frequency')
ylabel('Magnitude (dB)')
set(gca, 'Fontsize', 22)
title('PSD estimates and mean of noisy sinewave (dB)', 'Fontsize', 35)
figure
plot(normf, stdPSDs100dB, 'linewidth', 2);
xlabel('Normalised Frequency')
ylabel('Magnitude (dB)')
set(gca, 'Fontsize', 22)
title('STD of noisy sinewave PSD (dB)', 'Fontsize', 35)

% d: vary amount of points for periodogram resolution
M = 128;
numPoints = [30, 35, 50];

PSDd = zeros(3,M);

for i=1:length(numPoints)
    
    n = 0:numPoints(i);
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;
    x = [x zeros(1,M-length(x))];
    
    PSDd(i,:) = fft(x);
    
end

normf = linspace(0, 1, length(PSDd(1,:)));
figure
plot(normf,abs(PSDd(1,:)), 'linewidth', 2)
hold on
plot(normf,abs(PSDd(2,:)), 'linewidth', 2)
plot(normf,abs(PSDd(3,:)), 'linewidth', 2)
xlim([0 0.5])
xlabel('Normalised Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('n = 30','n = 35','n = 50')
title('Periodogram of 2 complex exponentials for varying n', 'Fontsize', 35)

% e MUSIC

M = 50;
N = 30;

n = 0:N;

x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n);
% x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+exp(1j*2*pi*0.38*n);

for i=1:M
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    xn = x + noise;
    
    [X,R] = corrmtx(xn,14,'mod');
    [S(:,i),F(:,i)] = pmusic(R,3,[ ],1,'corr');
end

meanS = mean(S, 2);

figure
plot(F,S,'c','linewidth',2);
hold on
plot(F(:,1), meanS, 'b', 'linewidth', 2);
set(gca,'xlim',[0.25 0.40]);
grid on; 
xlabel('Frequency (Hz)'); 
ylabel('Pseudospectrum');
set(gca, 'Fontsize', 22)
title('Periodogram estimates and mean of 2 complex exponentials', 'Fontsize', 35)

%% EX 1.4

N = 10000; %change to 10000 for 2nd part
n = 1:N;
sigma = 1;
x = randn(size(n));
% initialize first 4 elements
x(1) = x(1);
x(2) = 2.76*x(1) + x(2);
x(3) = 2.76*x(2) - 3.81*x(1) + x(3);
x(4) = 2.76*x(3) - 3.81*x(2) + 2.65*x(1) + x(4);

for i = 5:N
 x(i) = 2.76*x(i-1) - 3.81*x(i-2) + 2.65*x(i-3) - 0.92*x(i-4) + x(i);
end

x = x(5001:end);

p=2:14;
h = zeros(length(p), N/2);
aOrig = [2.76 -3.81 2.65 -0.92];
hOrig = freqz(1 , [1 -aOrig], N/2);
for i=1:length(p)
    [a,e] = aryule(x,p(i));
    [h(i,:),~] = freqz(e^(1/2),a,N/2);
end

normf = linspace(0,1,length(hOrig));
figure
plot(normf,pow2db(abs(h(2,:)').^2),'c', 'linewidth',2)
hold on
plot(normf,pow2db(abs(h(4,:)').^2),'b', 'linewidth',2)
hold on
plot(normf,pow2db(abs(h(8,:)').^2),'g', 'linewidth',2)
hold on
plot(normf,pow2db(abs(h(11,:)').^2),'y', 'linewidth',2)
hold on
plot(normf,pow2db(abs(hOrig).^2), 'k', 'linewidth',2) 
xlim([0 0.6])
xlabel('Normalised Frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('AR(2)','AR(4)','AR(8)','AR(11)','Original')
title('Spectrum estimation for varying orders', 'Fontsize', 35)

%% EX 1.5

load('RRI1.mat');
load('RRI2.mat');
load('RRI3.mat');
xRRI1 = zscore(xRRI1);
xRRI2 = zscore(xRRI2);
xRRI3 = zscore(xRRI3);
pdg1 = (1/length(xRRI1))*(abs(fft(xRRI1))).^2;

figure
AX = linspace(0, 4, length(pdg1)-1);
plot(AX,pow2db(pdg1(2:end)), 'linewidth', 2);
xlim([0 2])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Periodogram for Trial 1', 'Fontsize', 35)

pdg2 = (1/length(xRRI2))*(abs(fft(xRRI2))).^2;
figure
AX = linspace(0, 4, length(pdg2)-1);
plot(AX,pow2db(pdg2(2:end)), 'linewidth', 2);
xlim([0 2])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Periodogram for Trial 2', 'Fontsize', 35)

pdg3 = (1/length(xRRI3))*(abs(fft(xRRI3))).^2;
figure
AX = linspace(0, 4, length(pdg3)-1);
plot(AX,pow2db(pdg3(2:end)), 'linewidth', 2);
xlim([0 2])
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Periodogram for Trial 3', 'Fontsize', 35)

% No considerable improvement in increasing padding
xRRI1pad = [xRRI1, zeros(1,2000-length(xRRI1))];
xRRI2pad = [xRRI2, zeros(1,2000-length(xRRI2))];
xRRI3pad = [xRRI3, zeros(1,2000-length(xRRI3))];
%RRI length now 1000 samples
%xRRI = [xRRI1pad; xRRI2pad; xRRI3pad];
wlen = [50 100 200];
L=2048;
fs=fsRRI1;
[PSD1L1, f11] = pwelch(xRRI1pad, hamming(wlen(1)), 0, L, fs, 'onesided');
[PSD1L2, f12] = pwelch(xRRI1pad, hamming(wlen(2)), 0, L, fs, 'onesided');
[PSD1L3, f13] = pwelch(xRRI1pad, hamming(wlen(3)), 0, L, fs, 'onesided');
[PSD2L1, f21] = pwelch(xRRI2pad, hamming(wlen(1)), 0, L, fs, 'onesided');
[PSD2L2, f22] = pwelch(xRRI2pad, hamming(wlen(2)), 0, L, fs, 'onesided');
[PSD2L3, f23] = pwelch(xRRI2pad, hamming(wlen(3)), 0, L, fs, 'onesided');
[PSD3L1, f31] = pwelch(xRRI3pad, hamming(wlen(1)), 0, L, fs, 'onesided');
[PSD3L2, f32] = pwelch(xRRI3pad, hamming(wlen(2)), 0, L, fs, 'onesided');
[PSD3L3, f33] = pwelch(xRRI3pad, hamming(wlen(3)), 0, L, fs, 'onesided');

figure
hold on
plot(f11, pow2db(PSD1L1), 'linewidth', 2);
plot(f12, pow2db(PSD1L2), 'linewidth', 2);
plot(f13, pow2db(PSD1L3), 'linewidth', 2);
xlim([0 2])
grid on; grid minor;
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('window length = 50','window length = 100','window length = 200')
title('Averaged Periodogram for Trial 1', 'Fontsize', 35)
figure
hold on
plot(f21, pow2db(PSD2L1), 'linewidth', 2);
plot(f22, pow2db(PSD2L2), 'linewidth', 2);
plot(f23, pow2db(PSD2L3), 'linewidth', 2);
xlim([0 2])
grid on; grid minor;
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('window length = 50','window length = 100','window length = 200')
title('Averaged Periodogram for Trial 2', 'Fontsize', 35)
figure
hold on
plot(f31, pow2db(PSD3L1), 'linewidth', 2);
plot(f32, pow2db(PSD3L2), 'linewidth', 2);
plot(f33, pow2db(PSD3L3), 'linewidth', 2);
xlim([0 2])
grid on; grid minor;
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('window length = 50','window length = 100','window length = 200')
title('Averaged Periodogram for Trial 3', 'Fontsize', 35)

% c

PSD1 = pow2db(periodogram(xRRI1pad, hamming(length(xRRI1pad)), L, 'centered'));
PSD2 = pow2db(periodogram(xRRI2pad, hamming(length(xRRI2pad)), L, 'centered'));
PSD3 = pow2db(periodogram(xRRI3pad, hamming(length(xRRI3pad)), L, 'centered'));


p = [1 4 8 15 22 35];

for i = 1:length(p)   
    % Straight forward implementation    
    [a, e] = aryule(xRRI1pad, p(i));
    [A, w1] = freqz(e.^(1/2), a, length(xRRI1pad), fs);
    PSD1AR(i, :) = abs(A).^2;     
end

for i = 1:length(p)   
    % Straight forward implementation    
    [a, e] = aryule(xRRI2pad, p(i));
    [A, w2] = freqz(e.^(1/2), a, length(xRRI2pad), fs);
    PSD2AR(i, :) = abs(A).^2;     
end

for i = 1:length(p)   
    % Straight forward implementation    
    [a, e] = aryule(xRRI3pad, p(i));
    [A, w3] = freqz(e.^(1/2), a, length(xRRI3pad), fs);
    PSD3AR(i, :) = abs(A).^2;     
end

AX = linspace(0, 4, length(pdg1)-1);
f = (fs/2)*linspace(-1, 1, L);
figure
hold on
plot(f,PSD1, 'linewidth', 2);
plot(w1,pow2db(PSD1AR'), 'linewidth', 2);
xlim([0 2])
grid on; grid minor;
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Original','p=1','p=4','p=8','p=15','p=22','p=35')
title('AR spectrum estimation Trial 1', 'Fontsize', 35)

AX = linspace(0, 4, length(pdg2)-1);
figure
hold on
plot(f,PSD2, 'linewidth', 2);
plot(w1,pow2db(PSD2AR'), 'linewidth', 2);
xlim([0 2])
grid on; grid minor;
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Original','p=1','p=4','p=8','p=15','p=22','p=35')
title('AR spectrum estimation Trial 2', 'Fontsize', 35)

AX = linspace(0, 4, length(pdg3)-1);
figure
hold on
plot(f,PSD3, 'linewidth', 2);
plot(w1,pow2db(PSD3AR'), 'linewidth', 2);
xlim([0 2])
grid on; grid minor;
xlabel('Normalized frequency')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
legend('Original','p=1','p=4','p=8','p=15','p=22','p=35')
title('AR spectrum estimation Trial 3', 'Fontsize', 35)
%% EX 1.6

load 'PCAPCR.mat'

svdX = svd(X);
svdXN = svd(Xnoise);

figure
stem(svdX, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Singular Value Index')
ylabel('Singular Value Magnitude')
set(gca, 'Fontsize', 22)
title('Singular Values for signal X', 'Fontsize', 35)
figure
stem(svdXN, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Singular Value Index')
ylabel('Singular Value Magnitude')
set(gca, 'Fontsize', 22)
title('Singular Values for Xnoise', 'Fontsize', 35)

figure
stem((svdX-svdXN).^2, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Singular Value Index')
ylabel('Magnitude')
set(gca, 'Fontsize', 22)
title('Squared Error', 'Fontsize', 35)

% b perform PCA and compare difference

for rank=1:length(svdX)
    [U,S,V] = svd(Xnoise);
    
    U = U(:,1:rank);
    S = S(1:rank,1:rank);
    V = V(:,1:rank);

    Xrecon = U*S*V';

    e1(rank) = sum(sum((X-Xrecon).^2));
    e2(rank) = sum(sum((Xrecon-Xnoise).^2));
end

figure
plot(e1, 'linewidth', 2)
hold on
plot(e2, 'linewidth', 2)
xlabel('Rank')
ylabel('Error')
set(gca, 'Fontsize', 22)
legend('X - Xrecon','Xrecon - Xnoise')
title('Squared Error vs reconstruction rank', 'Fontsize', 35)

% c

Bols = (Xnoise'*Xnoise)\Xnoise'*Y;
Yols = Xnoise*Bols;

rank = 3;
[U,S,V] = svd(Xnoise);
U = U(:,1:rank);
S = S(1:rank,1:rank);
V = V(:,1:rank);
BPCR = V/S*U'*Y;
Ypcr = Xnoise*BPCR;

eOLS = sum(sum((Y-Yols).^2))/numel(Y);
ePCR = sum(sum((Y-Ypcr).^2))/numel(Y);

YolsT = Xtest*Bols;
YpcrT = Xtest*BPCR;

eOLST = sum(sum((Ytest-YolsT).^2))/numel(Ytest);
ePCRT = sum(sum((Ytest-YpcrT).^2))/numel(Ytest);

% d
addpath('./')
[YhatO, YO] = regval(Bols);
[YhatP, YP] = regval(BPCR);

eO = sum(sum((YO-YhatO).^2))/numel(YO);
eP = sum(sum((YP-YhatP).^2))/numel(YP);

