function Rho_k = Update_mu(M,N,K,h_k,g_k,H,Theta,eps_k,f_k)
    Rho_k=zeros(K,1);
    for k=1:K
        h_k_temp=reshape(h_k(k,:),N,1);
        g_k_temp=reshape(g_k(k,:),M,1);
        temp1=h_k_temp+H'*Theta(:,:,k)*g_k_temp;
        
        temp2=reshape(f_k(k,:),N,1);
        temp=real(eps_k(k)'*temp1'*temp2);
        Rho_k(k)=(temp^2+temp*sqrt(temp^2+4))/2;
    end
    
end

