function [R] = Active_Precoding(N,K,M,PRIS,sigma1,sigma2,Theta,f_k,h_k,g_k,H,beta_max)
iteration = 40;
obj = zeros(2*iteration, 1);
nu = ones(K, 1);
for i = 1:iteration
    mu = Update_mu(M,N,K,h_k,g_k,H,Theta,nu,f_k); 
    nu = MF_Active_Update_nu(K,N,M,mu,h_k,g_k,H,Theta,f_k,sigma1,sigma2); 
    [~,obj(2*i-1)] = MF_Active_Update_SINR(K,N,M,h_k,g_k,H,Theta,f_k,sigma1,sigma2);
    [v,Q] = MF_Active_Update_Q_v(K,N,M,mu,nu,h_k,g_k,H,f_k,sigma2);
    theta_temp1 = Active_Update_theta(M,K,N,v,Q,f_k,H,PRIS,sigma2,beta_max); 
    theta_temp2 = zeros(K,M);
    theta_phase = zeros(K,M);
    theta_amp = zeros(K,M);
    Theta = zeros(M,M,K);
    for k=1:K
        theta_temp2(k,:)=reshape(theta_temp1(M*(k-1)+1:1:M*k),M,1); 
        theta_phase(k,:) = angle(theta_temp2(k,:));
        theta_amp(k,:) = abs(theta_temp2(k,:));
        Theta(:,:,k) = diag(theta_amp(k,:).*exp(1j*theta_phase(k,:)));
    end

    [R,obj(2*i)] = MF_Active_Update_SINR(K,N,M,h_k,g_k,H,Theta,f_k,sigma1,sigma2); 
    
    if i>1
        if (obj(2*i)-obj(2*i-1))/obj(2*i-1)<0.01
            break;
        end
    end
end
end