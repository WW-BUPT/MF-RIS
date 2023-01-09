function [h] = Generate_Channel_h(M,N,dis)
h = zeros(M,N);
h = raylrnd(1,M,N);
for aa=1:M
    for bb=1:N
       h(aa,bb) = h(aa,bb)*exp(1j*2*pi*rand());
    end
end
h = sqrt(10^(-3)*dis^(-3.5))*h;
end
