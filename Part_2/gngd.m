function [xhat, error, A] = gngd(x, in, mu, order, eps, rho)
    N = length(x);
    A = zeros(order,N);
    error = zeros(1,N);
    for i=order+1:N
        inpast = in(i:-1:i-order+1);
        inpastprev = in(i-1:-1:i-order);
        xhat(i) = A(:,i)'*inpast';
        error(i) = x(i)-xhat(i);
        A(:,i+1) = A(:,i) + (1/(eps(i) + inpast*inpast'))*error(i)*inpast';
        
        top = eps(i)*eps(i-1)*inpast*inpastprev';
        bottom = (eps(i-1) + (inpastprev*inpastprev')).^2;
        eps(i+1) = eps(i) - rho*mu*(top/bottom);
    end
    A = A(:,2:end);
end