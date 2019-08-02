%% Ex 4

load('time-series.mat')
yOrig = y;
meanY = mean(y);
y = y-meanY; % remove mean from y

figure
plot(yOrig,'linewidth',2)
xlabel('n')
ylabel('Amplitude')
set(gca, 'Fontsize', 22)
title('Original time series signal', 'Fontsize', 35)

mu = 0.00001;
order = 4;

[yhat, error, A] = lms(y, mu, order);
MSE = mean(error.^2);

errorVar = var(error);
yhatVar = var(yhat);
predGain = 10*log10(yhatVar/errorVar);

figure
hold on
plot(y,'linewidth',2)
plot(yhat,'linewidth',2)
xlabel('n')
ylabel('Amplitude')
set(gca, 'Fontsize', 22)
legend('Zero mean Original','LMS prediction')
title('Time series signals', 'Fontsize', 35)

%--------  b  ---------

mu = 0.00001;
order = 4;

[yhat, error, A] = lmsD(y, mu, order,1);
MSE = mean(error.^2);

errorVar = var(error);
yhatVar = var(yhat);
predGain = 10*log10(yhatVar/errorVar);

figure
hold on
plot(y,'linewidth',2)
plot(yhat,'linewidth',2)
xlabel('n')
ylabel('Amplitude')
set(gca, 'Fontsize', 22)
legend(' Zero mean Original','LMS with tanh prediction')
title('Time series signals', 'Fontsize', 35)

%--------  c  ---------

mu = 0.00001;
order = 4;

[yhat, error, A] = lmsD(y, mu, order,40);
MSE = mean(error.^2);

errorVar = var(error);
yhatVar = var(yhat);
predGain = 10*log10(yhatVar/errorVar);

figure
hold on
plot(y,'linewidth',2)
plot(yhat,'linewidth',2)
xlabel('n')
ylabel('Amplitude')
set(gca, 'Fontsize', 22)
legend(' Zero mean Original','LMS with a.tanh prediction')
title('Time series signals', 'Fontsize', 35)

%--------  d  ---------

mu = 0.00001;
order = 4;

yAug = [1; yOrig];

[yhat, error, A] = lmsD(yAug, mu, order,50);
MSE = mean(error.^2);

errorVar = var(error);
yhatVar = var(yhat);
predGain = 10*log10(yhatVar/errorVar);

figure
hold on
plot(yOrig,'linewidth',2)
plot(yhat(2:end),'linewidth',2)
xlim([0 1000])
xlabel('n')
ylabel('Amplitude')
set(gca, 'Fontsize', 22)
legend('Original','LMS with biased a.tanh prediction')
title('Time series signals', 'Fontsize', 35)

%-------------  e  --------------

mu = 0.00001;
order = 4;

yAug = [1; yOrig];
Y = yAug(1:21,1);

A = zeros(order,length(Y));

for i=1:100
    [yhat, error, A] = lmsW(Y, mu, order,50, A);
    a = A(:,end);
    A = [zeros(4,length(Y))];
    A(:,order+1) = a;
end

yAug = [1; yOrig];

[yhat, error, A] = lmsW(yAug, mu, order,50, A);
MSE = mean(error.^2);

errorVar = var(error);
yhatVar = var(yhat);
predGain = 10*log10(yhatVar/errorVar);

figure
hold on
plot(yOrig,'linewidth',2)
plot(yhat(2:end),'linewidth',2)
xlabel('n')
ylabel('Amplitude')
set(gca, 'Fontsize', 22)
legend('Original','LMS with training prediction')
title('Time series signals', 'Fontsize', 35)
