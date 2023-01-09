function [R,obj] = MF_Active_Update_SINR(K,N,M,h_k,g_k,H,Theta,f_k,sigma1,sigma2)
gamma_k=zeros(K,1);

for k=1:K
    h_k_temp=reshape(h_k(k,:),N,1);
    g_k_temp=reshape(g_k(k,:),M,1);
    temp1=h_k_temp+H'*Theta(:,:,k)*g_k_temp;
    
    temp2=reshape(f_k(k,:),N,1);
    temp3=temp1'*temp2;
    
    temp4=0;
    for j=1:K
        f_j=reshape(f_k(j,:),N,1);
        temp4=temp4+norm(temp1'*f_j)^2;
    end
    temp4=temp4-norm(temp3)^2+norm(g_k_temp'*Theta(:,:,k)')^2*sigma2+sigma1;    
    gamma_k(k)=norm(temp3)^2/temp4;
end
R=0;
for k=1:K
    R=R+log2(1+gamma_k(k));
end
obj=R;
end