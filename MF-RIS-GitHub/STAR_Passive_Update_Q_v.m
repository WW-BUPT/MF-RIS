function [nu,Lam] = STAR_Passive_Update_Q_v(K,N,M,Rho_k,eps_k,h_k,g_k,H,f_k)

nu_k=zeros(K,M);  
Lam_k=zeros(K,M,M);

for k=1:K 
    temp_g_k=reshape(g_k(k,:),M,1);
    temp_h_k=reshape(h_k(k,:),N,1);
    temp_f_k=reshape(f_k(k,:),N,1);
    
    temp1=2*sqrt(1+Rho_k(k))*eps_k(k)'*diag(temp_g_k')*H*temp_f_k;
    
    temp2=zeros(N,N);
    for j=1:K
        temp3=reshape(f_k(j,:),N,1); 
        temp3=temp3*temp3';
        temp2=temp2+temp3;   
    end
    
    nu_k(k,:)=temp1-abs(eps_k(k))^2*diag(temp_g_k')*H*temp2*temp_h_k;
    
    Lam_k(k,:,:)=abs(eps_k(k))^2*diag(temp_g_k')*H*temp2*H'*diag(temp_g_k); 
end

nu=zeros(K*M,1);
Lam=zeros(K*M,K*M);

for n=1:M
    nu(n,1)=nu_k(1,n);
end

for n=(M+1):(2*M)
    nu(n,1)=nu_k(2,n-M);
end

bar_I = zeros(2,2,2);
bar_I(:,:,1)=diag([1,0]);
bar_I(:,:,2)=diag([0,1]);

for k=1:K
    Lam_temp=reshape(Lam_k(k,:,:),M,M);
    Lam=Lam+kron(bar_I(:,:,k),Lam_temp);
end

end

