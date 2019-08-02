%% EX 3.1
% signal length and ensemble size
N = 1000;
realisations = 100;
% params
b = [1.5 + 1i, 2.5 - 0.5i];
a = 1;
order = length(b);
mu = 0.1;
var = 1;

% variables
h = complex(zeros(realisations, order, N));
e = complex(zeros(realisations, N));
hAug = complex(zeros(realisations, order, N));
gAug = complex(zeros(realisations, order, N));
eAug = complex(zeros(realisations, N));

for i=1:realisations
   
    x = wgn(1,N,pow2db(var),'complex'); % noise
    y = compFilt(b,a,x); % signal
    
    [h(i,:,:), e(i,:)] = clms(y, x, mu, order);
    [hAug(i,:,:), gAug(i,:,:), eAug(i,:)] = aclms(y, x, mu, order);
    
end

figure
hold on
scatter(real(y),imag(y),'.','linewidth',2)
scatter(real(x),imag(x),'.','linewidth',2)
xlabel('Real')
ylabel('imag')
set(gca, 'Fontsize', 22)
legend('WLMA','WGN')
title('Distribution of x and y', 'Fontsize', 35)
xlim([-10 10])

e = abs(e).^2;
eAug = abs(eAug).^2;
e = squeeze(mean(e, 1));
eAug = squeeze(mean(eAug, 1));

figure
hold on
plot(pow2db(e),'linewidth',2)
plot(pow2db(eAug),'linewidth',2)
xlabel('n')
ylabel('Error (dB)')
set(gca, 'Fontsize', 22)
legend('CLMS','ACLMS')
title('Learning Curves for x and y', 'Fontsize', 35)
grid on
grid minor


%------------  b  --------------

load('./low-wind.mat')
vLow = v_east + 1i*v_north;
load('./medium-wind.mat')
vMed = v_east + 1i*v_north;
load('./high-wind.mat')
vHigh = v_east + 1i*v_north;

% Circularity
rhoLow = abs(mean((vLow).^2)/mean(abs(vLow).^2));
rhoMed = abs(mean((vMed).^2)/mean(abs(vMed).^2));
rhoHigh = abs(mean((vHigh).^2)/mean(abs(vHigh).^2));

% Center points
vLowR = mean(real(vLow));
vLowI = mean(imag(vLow));
vMedR = mean(real(vMed));
vMedI = mean(imag(vMed));
vHighR = mean(real(vHigh));
vHighI = mean(imag(vHigh));

figure
hold on
scatter(real(vLow),imag(vLow));
scatter(0,0, 100,'x','k','Linewidth',2)
scatter(vLowR,vLowI, 100,'x','b','Linewidth',2)
xlim([-0.4 0.4])
ylim([-0.4 0.4])
xlabel('Real')
ylabel('Imag')
set(gca, 'Fontsize', 22)
title(['Low Wind, rho=',num2str(rhoLow)])
grid on

figure
hold on
scatter(real(vMed),imag(vMed));
scatter(0,0, 100,'x','k','Linewidth',2)
scatter(vMedR,vMedI, 100,'x','b','Linewidth',2)
xlim([-2 2])
ylim([-2 2])
xlabel('Real')
ylabel('Imag')
set(gca, 'Fontsize', 22)
title(['Medium Wind, rho=',num2str(rhoMed)])
grid on

figure
hold on
scatter(real(vHigh),imag(vHigh));
scatter(0,0, 100,'x','k','Linewidth',2)
scatter(vHighR,vHighI, 100,'x','b','Linewidth',2)
xlim([-4 4])
ylim([-4 4])
xlabel('Real')
ylabel('Imag')
set(gca, 'Fontsize', 22)
title(['High Wind, rho=',num2str(rhoHigh)])
grid on

N = length(vLow);
mu = 0.001;
orderRange = 30;

% vLow
e = complex(zeros(orderRange,N));
eAug = complex(zeros(orderRange,N));

for r=1:orderRange
    y = vLow;
    x = [0; vLow(1:end-1)]';
    
    [h, e(r,:)] = clms(y, x, mu, r);
    [hAug, gAug, eAug(r,:)] = aclms(y, x, mu, r);
    
end

e = abs(e).^2;
eAug = abs(eAug).^2;
e = squeeze(mean(e, 2));
eAug = squeeze(mean(eAug, 2));

figure 
hold on
plot(pow2db(e),'Linewidth',2)
plot(pow2db(eAug),'Linewidth',2)
xlabel('Model Order')
ylabel('Error(dB)')
set(gca, 'Fontsize', 22)
legend('CLMS','ACLMS')
title(['MSPE for low wind'])
grid on


% vMed
e = complex(zeros(orderRange,N));
eAug = complex(zeros(orderRange,N));

for r=1:orderRange
    y = vMed;
    x = [0; vMed(1:end-1)]';
    
    [h, e(r,:)] = clms(y, x, mu, r);
    [hAug, gAug, eAug(r,:)] = aclms(y, x, mu, r);
    
end

e = abs(e).^2;
eAug = abs(eAug).^2;
e = squeeze(mean(e, 2));
eAug = squeeze(mean(eAug, 2));

figure 
hold on
plot(pow2db(e),'Linewidth',2)
plot(pow2db(eAug),'Linewidth',2)
xlabel('Model Order')
ylabel('Error(dB)')
set(gca, 'Fontsize', 22)
legend('CLMS','ACLMS')
title(['MSPE for medium wind'])
grid on


% vHigh
e = complex(zeros(orderRange,N));
eAug = complex(zeros(orderRange,N));

for r=1:orderRange
    y = vHigh;
    x = [0; vHigh(1:end-1)]';
    
    [h, e(r,:)] = clms(y, x, mu, r);
    [hAug, gAug, eAug(r,:)] = aclms(y, x, mu, r);
    
end

e = abs(e).^2;
eAug = abs(eAug).^2;
e = squeeze(mean(e, 2));
eAug = squeeze(mean(eAug, 2));

figure 
hold on
plot(pow2db(e),'Linewidth',2)
plot(pow2db(eAug),'Linewidth',2)
xlabel('Model Order')
ylabel('Error(dB)')
set(gca, 'Fontsize', 22)
legend('CLMS','ACLMS')
title(['MSPE for high wind'])
grid on



% ---------  c  -----------

% params
N = 1000;
f = 50;
fs = 10000;
Vmag = 1;
phi = 0;

% balanced
V = Vmag*ones(1,3);
Delta = zeros(1,3);

v = clarke(V, Delta, f, fs, phi, N);

figure
hold on
scatter(real(v),imag(v),'Linewidth',2);
xlim([-2 2])
ylim([-2 2])
xlabel('Real')
ylabel('Imag')
set(gca, 'Fontsize', 20)
title(['Balanced system'])
grid on

% unbalanced: magnitude

coeff = [1, 1, 1; 0.7, 1, 1.3; 0.4, 1, 1.6; 0.1, 1, 1.9];
v = zeros(4, N);
V = Vmag*ones(1,3);
Delta = zeros(1,3);

for i=1:size(coeff,1)
    V2 = V.*coeff(i,:);
    v(i,:) = clarke(V2, Delta, f, fs, phi, N);
end

figure
hold on
scatter(real(v(1,:)),imag(v(1,:)),'Linewidth',2)
scatter(real(v(2,:)),imag(v(2,:)),'Linewidth',2)
scatter(real(v(3,:)),imag(v(3,:)),'Linewidth',2)
scatter(real(v(4,:)),imag(v(4,:)),'Linewidth',2)
xlim([-2 2])
ylim([-2 2])
xlabel('Real')
ylabel('Imag')
set(gca, 'Fontsize', 20)
legend('V=[1,1,1]','V=[0.7,1,1.3]','V=[0.4,1,1.6]','V=[0.1,1,1.9]')
title(['Unbalanced system: Voltage'])
grid on

% unbalanced: phase

coeff = [0, 0, 0; 0, 0.2, 0.3; 0, -1/2, 1/2; 0, 0.3, -0.3; 0, -0.3, 0.3];
v = zeros(4, N);
V = Vmag*ones(1,3);
Delta = ones(1,3);

for i=1:size(coeff,1)
    delta = Delta.*coeff(i,:)*pi;
    v(i,:) = clarke(V, delta, f, fs, phi, N);
end

figure
hold on
scatter(real(v(1,:)'),imag(v(1,:)'),'Linewidth',2)
scatter(real(v(2,:)'),imag(v(2,:)'),'Linewidth',2)
scatter(real(v(3,:)'),imag(v(3,:)'),'Linewidth',2)
scatter(real(v(4,:)'),imag(v(4,:)'),'Linewidth',2)
xlim([-2 2])
ylim([-2 2])
xlabel('Real')
ylabel('Imag')
set(gca, 'Fontsize', 20)
legend('Delta=[0,0,0]','Delta=[0,0.2,0.3]','Delta=[0,-0.5,0.5]','Delta=[0,0.3,-0.3]')
title(['Unbalanced system: Phase'])
grid on



%--------------  e  -------------

% params
N = 1000;
f = 50;
fs = 5000;
Vmag = 1;
phi = 0;
mu = 0.01;
order = 1;

% balanced
V = Vmag*ones(1,3);
Delta = zeros(1,3);

v = clarke(V, Delta, f, fs, phi, N);
% Circularity
rho = abs(mean((v).^2)/mean(abs(v).^2));

y = v;
x = [0; v(1:end-1)'].';

[h, e] = clms(y, x, mu, order);
[hAug, gAug, eAug] = aclms(y, x, mu, order);

for i=1:N
    fHat(i) = (fs/(2*pi))*atan(sqrt(imag(h(i)).^2-abs(0).^2)/real(h(i)));
    fHatAug(i) = (fs/(2*pi))*atan(sqrt(imag(hAug(i)).^2-abs(gAug(i)).^2)/real(hAug(i)));
end

figure 
hold on
plot(fHat,'Linewidth',2)
plot(abs(fHatAug),'Linewidth',2)
xlabel('n')
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 20)
legend('CLMS','ACLMS')
title(['Balanced System'])
grid on

% unbalanced: magnitude
V = Vmag*ones(1,3).*[0.6, 1, 1.4];
Delta = zeros(1,3);

v = clarke(V, Delta, f, fs, phi, N);
% Circularity
rho = abs(mean((v).^2)/mean(abs(v).^2));

y = v;
x = [0; v(1:end-1)'].';

[h, e] = clms(y, x, mu, order);
[hAug, gAug, eAug] = aclms(y, x, mu, order);

for i=1:N
    fHat(i) = (fs/(2*pi))*atan(sqrt(imag(h(i)).^2-abs(0).^2)/real(h(i)));
    fHatAug(i) = (fs/(2*pi))*atan(sqrt(imag(hAug(i)).^2-abs(gAug(i)).^2)/real(hAug(i)));
end

figure 
hold on
plot(fHat,'Linewidth',2)
plot(abs(fHatAug),'Linewidth',2)
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 20)
legend('CLMS','ACLMS')
title(['Unbalanced System: Voltage'])
grid on

% unbalanced: phase
V = Vmag*ones(1,3);
Delta = ones(1,3).*[0, 0.3, 0.2]*pi;

v = clarke(V, Delta, f, fs, phi, N);
% Circularity
rho = abs(mean((v).^2)/mean(abs(v).^2));

y = v;
x = [0; v(1:end-1)'].';

[h, e] = clms(y, x, mu, order);
[hAug, gAug, eAug] = aclms(y, x, mu, order);

for i=1:N
    fHat(i) = (fs/(2*pi))*atan(sqrt(imag(h(i)).^2-abs(0).^2)/real(h(i)));
    fHatAug(i) = (fs/(2*pi))*atan(sqrt(imag(hAug(i)).^2-abs(gAug(i)).^2)/real(hAug(i)));
end

figure 
hold on
plot(fHat,'Linewidth',2)
plot(abs(fHatAug),'Linewidth',2)
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 20)
legend('CLMS','ACLMS')
title(['Unbalanced System: Phase'])
grid on

%% Ex 3.2

% Params
N = 1500;
fs = 10000;
var = 0.05;
order = [1 5 10 20];

[phi, dphi] = phiSignal(N);
phase = ((2*pi)/fs).*phi;
eta = wgn(1,N,pow2db(var),'complex');
y = exp(1i*phase) + eta;

figure
plot(dphi,'Linewidth',2)
ylim([0 500])
xlabel('n')
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 20)
title(['Frequency of FM signal'])
grid on

a = zeros(N, length(order));
w = zeros(N, length(order));
for p=1:length(order)
    coeffs = aryule(y, order(p));
    [a(:,p), w(:,p)] = freqz(1, coeffs, N, fs);
end

figure
hold on
plot(w, mag2db(abs(a)),'Linewidth',2)
xlim([0 1000])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
set(gca, 'Fontsize', 20)
legend('AR(1)','AR(5)','AR(10)','AR(20)')
title(['Power Spectrum of FM signal'])
grid on

p = [1 10 20];
const1 = aryule(y(1:500), p(1));
const10 = aryule(y(1:500), p(2));
const20 = aryule(y(1:500), p(3));
lin1 = aryule(y(501:1000), p(1));
lin10 = aryule(y(501:1000), p(2));
lin20 = aryule(y(501:1000), p(3));
quad1 = aryule(y(1001:1500), p(1));
quad10 = aryule(y(1001:1500), p(2));
quad20 = aryule(y(1001:1500), p(3));

[hConst1, wConst1] = freqz(1, const1, N, fs);
[hConst10, wConst10] = freqz(1, const10, N, fs);
[hConst20, wConst20] = freqz(1, const20, N, fs);
[hLin1, wLin1] = freqz(1, lin1, N, fs);
[hLin10, wLin10] = freqz(1, lin10, N, fs);
[hLin20, wLin20] = freqz(1, lin20, N, fs);
[hQuad1, wQuad1] = freqz(1, quad1, N, fs);
[hQuad10, wQuad10] = freqz(1, quad10, N, fs);
[hQuad20, wQuad20] = freqz(1, quad20, N, fs);

figure
hold on
plot(wConst1, hConst1,'Linewidth',2)
plot(wConst10, hConst10,'Linewidth',2)
plot(wConst20, hConst20,'Linewidth',2)
xlim([0 500])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
set(gca, 'Fontsize', 20)
legend('AR(1)','AR(10)','AR(20)')
title(['Power Spectrum of constant section'])
grid on

figure
hold on
plot(wLin1, hLin1,'Linewidth',2)
plot(wLin10, hLin10,'Linewidth',2)
plot(wLin20, hLin20,'Linewidth',2)
xlim([0 500])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
set(gca, 'Fontsize', 20)
legend('AR(1)','AR(10)','AR(20)')
title(['Power Spectrum of linear section'])
grid on

figure
hold on
plot(wQuad1, hQuad1,'Linewidth',2)
plot(wQuad10, hQuad10,'Linewidth',2)
plot(wQuad20, hQuad20,'Linewidth',2)
xlim([0 600])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
set(gca, 'Fontsize', 20)
legend('AR(1)','AR(10)','AR(20)')
title(['Power Spectrum of exponential section'])
grid on


%-------------  b  -------------

% Params
L = 1024;
N = 1500;
fs = 2000;
var = 0.05;
mu = 0.3;

[phi, dphi] = phiSignal(N);
phase = ((2*pi)/fs).*phi;
eta = wgn(1,N,pow2db(var),'complex');
y = exp(1i*phase) + eta;

%clms
a = zeros(1, N);
H = zeros(L, N-1);
order=1;

x = [0; y(1:end-1)'].'; 
[a, error] = clms(y, x, mu, order);
for n=1:N
    [h, w] = freqz(1, [1; -conj(a(n))], L);
    H(:, n) = abs(h).^2;
end
% Remove outliers in H
medianH = 50*median(median(H));
H(H > medianH) = medianH;

figure
surf(1:N, (w/pi)*1000, H, 'LineStyle','none');
view(2)
xlabel('n')
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 20)
title(['FM signal estimated spectrum ',num2str(mu)])
grid on


%% Ex 3.3

% Parameters
N = 1500;
K = 500;
fs = 2000;
var = 0.05;
mu = 1;
gamma = 0.2; % change for leaky clms

[phi, dphi] = phiSignal(N);
phase = ((2*pi)/fs).*phi;
eta = wgn(1,N,pow2db(var),'complex');
y = exp(1i*phase) + eta;

w = zeros(K, N);
e = zeros(1, N);
input = (1/K)*exp(1i*2*pi*(0:N-1)'*(0:K-1)/K);

[w, e] = clmsLeaky(y, input, mu, gamma);

figure
f = (K-1:-1:0)*fs/K;
surf(1:N, f, abs(w), 'LineStyle','none');
view(2)
ylim([0 600])
colorbar
xlabel('n')
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 20)
title(['FM signal estimated spectrum gamma=',num2str(gamma)])% gamma=',num2str(gamma)])

%--------------  d  ----------------
%%
load('./EEG_Data/EEG_Data_Assignment1.mat')

% Params
N = 1200;
K = 16384;
mu = 1;
gamma = 0;
start = 36000;

y = POz(start:start+N-1);

w = zeros(K, N);
e = zeros(1, N);
input = (1/K)*exp(1i*2*pi*(0:N-1)'*(0:K-1)/K);

[w, e] = clmsLeaky(y, input, mu, gamma);

figure
f = (K-1:-1:0)*fs/K;
surf(1:N, f, abs(w), 'LineStyle','none');
view(2)
ylim([0 60])
xlim([0 1000])
colorbar
xlabel('n')
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 20)
title(['EEG POz signal estimated spectrum'])
