function y = compFilt(b,a,x)
    
    N = length(x);
    y = complex(zeros(1,N));
    
    for i=2:N
        y(i) = a*x(i) + b(1)*x(i-1) + b(2)*conj(x(i-1));
    end

end