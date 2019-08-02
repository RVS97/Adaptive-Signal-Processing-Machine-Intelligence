function [h, error] = clms(x, in, mu, order)
    N = length(x);
    h = complex(zeros(order,N));%+1;
    error = complex(zeros(1,N));
    for i=order+1:N
        inpast = in(i:-1:i-order+1);
        xhat(i) = h(:,i)'*inpast';
        error(i) = x(i)-xhat(i);
        h(:,i+1) = h(:,i) + mu*conj(error(i))*inpast';
    end
    h = h(:,2:end);
    
end