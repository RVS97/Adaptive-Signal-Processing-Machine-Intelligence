function [h, error] = clmsLeaky(x, in, mu, gamma)
    N = length(x);
    h = zeros(size(in'));
    error = zeros(1,N);
    for i=1:N
        input = in(i,:);
        xhat(i) = h(:,i)'*input';
        error(i) = x(i)-xhat(i);
        h(:,i+1) = (1-gamma*mu)*h(:,i) + mu*conj(error(i))*input';
    end
    h = h(:,2:end);
end