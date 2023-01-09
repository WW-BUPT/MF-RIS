function [g] = Generate_Channel_g(M,N,dis)
g = zeros(M,N);
g = raylrnd(1,M,N);
for aa=1:M
    for bb=1:N
       g(aa,bb) = g(aa,bb)*exp(1j*2*pi*rand());
    end
end
g = sqrt(10^(-3)*dis^(-3))*g;
end

