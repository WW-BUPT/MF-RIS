function eps_k=STAR_Passive_Update_nu(K,N,M,Rho_k,h_k,g_k,H,Theta,f_k,sigma1)
eps_k=zeros(K,1);
for k=1:K
    h_k_temp=reshape(h_k(k,:),N,1);
    g_k_temp=reshape(g_k(k,:),M,1);
    temp1=h_k_temp+H'*Theta(:,:,k)*g_k_temp;
    
    temp2=reshape(f_k(k,:),N,1);
    temp3=sqrt(1+Rho_k(k))*temp1'*temp2;
    
    temp4=0;
    for j=1:K
        f_j=reshape(f_k(j,:),N,1);
        temp4=temp4+norm(temp1'*f_j)^2;
    end
    temp4=temp4+sigma1;    
    eps_k(k)=temp3/temp4;
end


end

