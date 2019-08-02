function [xhat, error, A] = lmsW(x, mu, order, alpha, A)
    N = length(x);
    xhat = zeros(size(x));
    error = zeros(1,N);
    for i=order+1:N
        xpast = x(i-1:-1:i-order);
        xhat(i) = alpha*tanh(A(:,i)'*xpast);
        error(i) = x(i)-xhat(i);
        A(:,i+1) = A(:,i) + mu*error(i)*xpast;
    end
    A = A(:,2:end);
end