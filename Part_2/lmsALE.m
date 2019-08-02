function [xhat, error, A] = lmsALE(x, mu, M, delta)
    N = length(x);
    A = zeros(M,N);
    xhat = zeros(size(x));
    error = zeros(1,N);
    for i=delta+M:N
        xpast = x(i-delta:-1:i-delta-M+1);
        xhat(i) = A(:,i)'*xpast';
        error(i) = x(i)-xhat(i);
        A(:,i+1) = A(:,i) + mu*error(i)*xpast';
    end
    A = A(:,2:end);
end