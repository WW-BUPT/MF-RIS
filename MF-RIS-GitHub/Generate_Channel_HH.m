function [H] = Generate_Channel_HH(M,N,dis)
H = zeros(M,N);
H = raylrnd(1,M,N);
for aa=1:M
    for bb=1:N
       H(aa,bb) = H(aa,bb)*exp(1j*2*pi*rand());
    end
end
H = sqrt(10^(-3)*dis^(-2.8))*H;
end

