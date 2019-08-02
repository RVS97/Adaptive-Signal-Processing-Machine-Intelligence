function [xhat, error, A] = lmsG(x, in, mu, order)
    N = length(x);
    A = zeros(order,N);
    error = zeros(1,N);
    for i=order+1:N
        inpast = in(i:-1:i-order+1);
        xhat(i) = A(:,i)'*inpast';
        error(i) = x(i)-xhat(i);
        A(:,i+1) = A(:,i) + mu*error(i)*inpast';
    end
    A = A(:,2:end);
end