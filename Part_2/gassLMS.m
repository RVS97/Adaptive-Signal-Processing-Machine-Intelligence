function [xhat, error, A, MU] = gassLMS(x, in, mu, order, rho, type)
    alpha = 0.8;

    N = length(x);
    A = zeros(order,N);
    error = zeros(1,N);
    xhat = zeros(1,N);
    psi = zeros(order,N);
    MU = mu*ones(1,N);
    for i=order+1:N
        inpast = in(i:-1:i-order+1);
        xhat(i) = A(:,i)'*inpast';
        error(i) = x(i)-xhat(i);
        A(:,i+1) = A(:,i) + MU(i)*error(i)*inpast';
        MU(i+1) = MU(i) + rho*error(i)*inpast*psi(:,i);
        switch type
            case 'ben'
                M = MU(i)*inpast'*inpast;
                psi(:,i+1) = (eye(size(M))-M)*psi(:,i) + error(i)*inpast';
            case 'ang'
                psi(:,i+1) = alpha*psi(:,i) + error(i)*inpast';
            case 'mat'
                psi(:,i+1) = error(i)*inpast';
        end
    end
    A = A(:,2:end);
end