%% EX 2.1
% b
N=1000;
realisations = 1; % change for 1 or 100 realisations
a = [0.1, 0.8];
order = length(a);
var = 0.25;

% Generate Model
model = arima('Constant', 0, 'AR', a, 'Variance', var);
X = simulate(model, N, 'NumPaths', realisations);
%LMS
mu=0.05;
error = zeros(realisations,N);
A1 = zeros(length(a),N,realisations);
for j=1:realisations
    [xhat, error(j,:), A1(:,:,j)] = lms(X(:,j), mu, order);
    
end
error = error.^2;
error = mean(error,1);

mu=0.01;
error2 = zeros(realisations,N);
A2 = zeros(length(a),N,realisations);
for j=1:realisations
    [xhat, error2(j,:), A2(:,:,j)] = lms(X(:,j), mu, order);
end
error2 = error2.^2;
error2 = mean(error2,1);


figure
plot(pow2db(error),'linewidth',2)
hold on
plot(pow2db(error2),'linewidth',2)
xlabel('n')
ylabel('Error')
set(gca, 'Fontsize', 22)
legend('mu=0.05','mu=0.01')
title('Squared Error 100 real.', 'Fontsize', 35)

% c

e1 = mean(error(500:end));
e2 = mean(error2(500:end));
e1 = (e1-var)/var
e2 = (e2-var)/var

% d
figure
plot(mean(A1,3)','linewidth',2)
xlabel('n')
ylabel('Coeff')
set(gca, 'Fontsize', 22)
legend('a1','a2')
title('Coeff evolution mu=0.05, 1 real.', 'Fontsize', 35)
grid on

figure 
plot(mean(A2,3)','linewidth',2)
xlabel('n')
ylabel('Coeff')
set(gca, 'Fontsize', 22)
legend('a1','a2')
title('Coeff evolution mu=0.01, 1 real.', 'Fontsize', 35)
grid on
% f

%LMS
mu=0.05;
gamma = 0.5;
error = zeros(realisations,N);
A1 = zeros(length(a),N,realisations);
for j=1:realisations
    [xhat, error(j,:), A1(:,:,j)] = leakylms(X(:,j), mu, gamma, order);
    
end
mu=0.01;
gamma = 0.5;
error2 = zeros(realisations,N);
A2 = zeros(length(a),N,realisations);
for j=1:realisations
    [xhat, error2(j,:), A2(:,:,j)] = leakylms(X(:,j), mu, gamma, order);
end

figure
plot(mean(A1,3)','linewidth',2)
xlabel('n')
ylabel('Coeff')
set(gca, 'Fontsize', 22)
legend('a1','a2')
title('Coeff evol. mu=0.05 gamma=0.5, 100 real.', 'Fontsize', 35)
grid on
grid minor

figure 
plot(mean(A2,3)','linewidth',2)
xlabel('n')
ylabel('Coeff')
set(gca, 'Fontsize', 22)
legend('a1','a2')
title('Coeff evol. mu=0.01 gamma=0.5, 100 real.', 'Fontsize', 35)
grid on
grid minor

%% EX 2.2

N = 1000; % data length
realisations = 100; % number of realisations

% filter
b = [1 0.9];
a = 1;
order = length(b);

% params
rho = 0.005;
mu = [0.01, 0.1];

A = zeros(5, order, N, realisations);
error = zeros(5, N, realisations);
for i=1:realisations
    eta = normrnd(0, sqrt(0.5), 1, N);
    x = filter(b,a,eta);
    
    [xhat, error(1, :, i), A(1, :, :, i)] = lmsG(x, eta, mu(1), order);
    [xhat, error(2, :, i), A(2, :, :, i)] = lmsG(x, eta, mu(2), order);
    [xhat, error(3, :, i), A(3, :, :, i)] = gassLMS(x, eta, 0, order, rho, 'ben');
    [xhat, error(4, :, i), A(4, :, :, i)] = gassLMS(x, eta, 0, order, rho, 'ang');
    [xhat, error(5, :, i), A(5, :, :, i)] = gassLMS(x, eta, 0, order, rho, 'mat');
end

error = mean(error,3);
A = mean(A,4);

realA = b(2)*ones(1000,1);
wE1 = realA - squeeze(A(1,2,:));
wE2 = realA - squeeze(A(2,2,:));
wE3 = realA - squeeze(A(3,2,:));
wE4 = realA - squeeze(A(4,2,:));
wE5 = realA - squeeze(A(5,2,:));

figure
plot(wE1,'linewidth',2)
hold on
plot(wE2,'linewidth',2)
plot(wE3,'linewidth',2)
plot(wE4,'linewidth',2)
plot(wE5,'linewidth',2)
legend('LMS 0.01','LMS 0.1','ben','ang','mat')
xlabel('n')
ylabel('Error')
set(gca, 'Fontsize', 22)
title('Weight Error Curves', 'Fontsize', 35)

figure
plot(pow2db(wE1.^2),'linewidth',2)
hold on
plot(pow2db(wE2.^2),'linewidth',2)
plot(pow2db(wE3.^2),'linewidth',2)
plot(pow2db(wE4.^2),'linewidth',2)
plot(pow2db(wE5.^2),'linewidth',2)
legend('LMS 0.01','LMS 0.1','ben','ang','mat')
xlabel('n')
ylabel('Error')
set(gca, 'Fontsize', 22)
title('Squared Error Curves', 'Fontsize', 35)
% DO VARYING rho???


% c

N = 1000; % data length
realisations = 100; % number of realisations

% filter
b = [1 0.9];
a = 1;
order = length(b);

% params
rho = 0.001;
mu = 0.01;

A = zeros(2, order, N, realisations);
error = zeros(2, N, realisations);
eps = ones(1,N); 
for i=1:realisations
    eta = normrnd(0, sqrt(0.5), 1, N);
    x = filter(b,a,eta);
    
    [xhat, error(1, :, i), A(1, :, :, i)] = gassLMS(x, eta, mu, order, rho, 'ben');
    [xhat, error(2, :, i), A(2, :, :, i)] = gngd(x, eta, mu, order, eps, rho);
end

error = mean(error,3);
A = mean(A,4);

realA = b(2)*ones(1,1000);
wBenError = realA - squeeze(A(1,2,:));
wGngdError = realA - squeeze(A(2,2,:));

figure
plot(wBenError(:,1), 'b','linewidth',2)
hold on
plot(wGngdError(:,1), 'k','linewidth',2)
xlim([0 200])
legend('Ben','GNGD')
xlabel('n')
ylabel('Error')
set(gca, 'Fontsize', 22)
title('Squared Error Curves rho=0.001 mu=0.01', 'Fontsize', 35)
grid on
grid minor


%% EX 2.3

N = 1000; % data length
realisations = 100; % number of realisations

% filter
b = [1 0 0.5];
a = 1;
order = length(b);

% params
mu = 0.01;
delta = [1 2 3 4];
M = 5;

xhat = zeros(length(delta),N,realisations);
error = zeros(length(delta),N,realisations);
A = zeros(length(delta),M,N,realisations);
s = zeros(realisations,N);
for d=1:length(delta)
    
    x = sin(0.01*pi*(1:N));
    for i=1:realisations
        v = normrnd(0, sqrt(1), 1, N);
        eta = filter(b, a, v);
        s(i,:) = x + eta;
        
        [xhat(d,:,i), error(d, :, i), A(d, :, :, i)] = lmsALE(s(i,:), mu, M, delta(d));
    end
    MSPE(d) = mean(mean((repmat(x', 1, realisations)-squeeze(xhat(d,:,:))).^2));
end

figure
plot(s', 'r','linewidth',2);
hold on
plot(squeeze(xhat(3,:,:)), 'b','linewidth',2); % vary for different plots
plot(x', 'k','linewidth',2);
xlabel('n')
ylabel('Signal Amplitude')
set(gca, 'Fontsize', 22)
title('MSPE=0.3077', 'Fontsize', 35)

figure
plot(mean(xhat(3,:,:),3),'linewidth',2)
xlabel('n')
ylabel('Denoised signal')
set(gca, 'Fontsize', 22)
title('Ensemble average Delta=4', 'Fontsize', 35)

% b 

N = 1000; % data length
realisations = 100; % number of realisations

% filter
b = [1 0 0.5];
a = 1;
order = length(b);

% params
mu = 0.01;
delta = 1:25;
M = [5 10 15 20];


MSPE = zeros(length(M),length(delta));
for m=1:length(M)
    xhat = zeros(length(delta),N,realisations);
    error = zeros(length(delta),N,realisations);
    A = zeros(length(delta),M(m),N,realisations);
    s = zeros(realisations,N);
    for d=1:length(delta)

        x = sin(0.01*pi*(1:N));
        for i=1:realisations
            v = normrnd(0, sqrt(1), 1, N);
            eta = filter(b, a, v);
            s(i,:) = x + eta;

            [xhat(d,:,i), error(d, :, i), A(d, :, :, i)] = lmsALE(s(i,:), mu, M(m), delta(d));
        end
        MSPE(m,d) = mean(mean((repmat(x', 1, realisations)-squeeze(xhat(d,:,:))).^2));
    end
end

figure
plot(MSPE','linewidth',2)
xlabel('Delta')
ylabel('MSPE')
set(gca, 'Fontsize', 22)
legend('M=5','M=10','M=15','M=20')
title('MSPE vs delta and M', 'Fontsize', 35)

figure
plot(M,MSPE(:,3),'linewidth',2)
xlabel('order')
ylabel('MSPE')
set(gca, 'Fontsize', 22)
title('MSPE vs order', 'Fontsize', 35)

% c

b = [1 0 0.5];
a = 1;
order = length(b);

% params
mu = 0.01;
delta = [3];
M = 5;

xhat = zeros(N,realisations);
error = zeros(N,realisations);
A = zeros(M,N,realisations);
s = zeros(realisations,N);

    
x = sin(0.01*pi*(1:N));
for i=1:realisations
    v = normrnd(0, sqrt(1), 1, N);
    eta = filter(b, a, v);
    s(i,:) = x + eta;

    [xhat(:,i), error(:, i), A(:, :, i)] = lmsG(s(i,:), eta, mu, M); %x, in, mu, order)
end
MSPE = mean(mean((repmat(x', 1, realisations)-xhat(:,:)).^2));


figure
plot(s', 'r','linewidth',2);
hold on
plot(error(:,:), 'b','linewidth',2); % vary for different plots
plot(x', 'k','linewidth',2);
xlabel('n')
ylabel('Signal Amplitude')
set(gca, 'Fontsize', 22)
title('ANC Delta=3 M=5', 'Fontsize', 35)

%%
% d

load('.\EEG_Data\EEG_Data_Assignment2.mat')

x = detrend(POz);

N = length(x);
% signal generation
f = 50/fs;
var = 1;
sine = sin(2*pi*f*(1:N));
w = normrnd(0, sqrt(var), 1, N);
eps = sine + w;

% Spectogram
K = 2^12;
L = 16384;
overlap = 0.75;
figure
spectrogram(x, hamming(K), floor(overlap*K),L, fs,'yaxis');
ylim([0 100])
xlabel('Time (min)')
ylabel('Frequency (Hz)')
set(gca, 'Fontsize', 22)
title(['Original signal'])

% varying mu and M : comment/uncomment accordingly
% mu = [0.001, 0.005, 0.01, 0.03];
% M = 20;
mu = 0.01;
M = [5 10 20 30];

l = length(M); % change mu-M accordingly

% change init according to what tested
xhat = zeros(N,l);
error = zeros(N,l);
% PSD = zeros(L,l);
for i=1:l
    
%     A = zeros(M(i),N);
%     [xhat(:,i), error(:, i), A(:, :)] = lmsG(x, eps, mu, M(i)); %x, in, mu, order)

    A = zeros(M(i),N);
    [xhat(:,i), error(:, i), A(:, :)] = lmsG(x, eps, mu, M(i)); %x, in, mu, order)
    
    figure
    spectrogram(squeeze(error(:,i)), hamming(K), floor(overlap*K),L, fs,'yaxis');
    ylim([0 100])
    xlabel('Time (min)')
    ylabel('Frequency (Hz)')
    set(gca, 'Fontsize', 22)
    title(['M=', num2str(M(i)), ' \mu=', num2str(mu)])
    
    xHat = reshape(error(:,i)-mean(x), 5*fs, []);
    PSD(:,i) = mean(periodogram(xHat,hamming(length(xHat)),5*fs,fs), 2);
%     PSD(:,i) = mean(periodogram(xHat,hamming(length(xHat)), [], L, 'centered'), 2);   
end

xr = reshape(x, 5*fs, []);
% PSDOrig = mean(periodogram(xr, [], L, 'centered'), 2);
PSDOrig = mean(periodogram(xr,hamming(length(xHat)),5*fs,fs), 2);
fp = (fs/2)*linspace(0, 1, length(PSDOrig));

figure
plot(fp,10*log10(PSDOrig), 'linewidth', 2);
hold on
plot(fp,10*log10(PSD(:,1)), 'linewidth', 2);
plot(fp,10*log10(PSD(:,2)), 'linewidth', 2);
plot(fp,10*log10(PSD(:,3)), 'linewidth', 2);
plot(fp,10*log10(PSD(:,4)), 'linewidth', 2);
xlim([0 80])
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
set(gca, 'Fontsize', 22)
legend('original','M=5','M=10','M=20','M=30')
title(['Periodogram for varying M'])


