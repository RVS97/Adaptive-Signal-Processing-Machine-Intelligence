function [phi, dphi] = phiSignal(N)
    part1 = 100*ones(1,500);
    part2 = 100 + ((501:1000)-500)./2;
    part3 = 100 + (((1001:1500) - 1000)./25).^2;
    dphi = [part1, part2, part3];
    phi = cumsum(dphi);
end