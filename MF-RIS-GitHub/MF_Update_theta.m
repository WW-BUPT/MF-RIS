function theta = MF_Update_theta(M,K,N,nu,Lam,f_k,H,PRIS,sigma2,beta_max)

bar_I = zeros(2,2,2);
bar_I(:,:,1)=diag([1,0]);
bar_I(:,:,2)=diag([0,1]);

U=zeros(K*M,K*M);
for k=1:K
    f_k_temp=reshape(f_k(k,:),N,1);
    beta_k = H*f_k_temp;
    U_temp=diag(abs(beta_k).^2)+sigma2*eye(M);
    U=U+kron(bar_I(:,:,k),U_temp) ;
end

Lam=0.5*(Lam+Lam'); 
U=0.5*(U+U');

R_n = zeros(2*M,2*M,2*M);
for n_m=1:(2*M)
   for n_out=n_m
        for n_in=n_m
        R_n(n_in,n_in,n_m)=1;
        end
    end
end

cvx_begin quiet
    variable theta(K*M,1) complex;
    maximize(-(theta')*Lam*theta+2*real((theta')*nu)) 
    subject to
    theta'*U*theta <= PRIS;
    for n = 1: (2*M)
        theta'*R_n(:,:,n)*theta <= beta_max;
    end
cvx_end    

end