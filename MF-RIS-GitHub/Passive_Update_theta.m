function theta = Passive_Update_theta(M,nu,Lam)

Lam=0.5*(Lam+Lam');

R_n = zeros(2*M,2*M,2*M);
for n_m=1:(2*M)
   for n_out=n_m
        for n_in=n_m
        R_n(n_in,n_in,n_m)=1;
        end
    end
end

cvx_begin quiet
    variable theta(2*M,1) complex;
    maximize(-(theta')*Lam*theta+2*real((theta')*nu)) 
    subject to
    for n = 1: (2*M)
        theta'*R_n(:,:,n)*theta <= 1;
    end
    for n=1:(M/2)
        theta(n,1)==0;
    end
    for n=((3/2)*M+1):(2*M)
        theta(n,1)==0;
    end
cvx_end    

end