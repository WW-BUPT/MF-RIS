clear all;
User_Position = 20;
Iteration = 40;
N = 2; K = 2;  
M_ite = 10:4:50;
PBS = db2pow(9); PRIS = db2pow(9);
sigma1 = 1e-11;
sigma2 = 1e-11; 

beta_max = 1000; 

R3 = zeros(length(M_ite), Iteration); R2 = R3; R1 = R2;  R0 = R1;

tic;
F= exp(1j*2*pi*rand(K*N,1))*sqrt(PBS/K/N);
f_k=zeros(K,N);
for k=1:K
    f_k(k,:)=reshape(F(N*(k-1)+1:1:N*k),N,1);  
end

for u = 1:length(M_ite)
      parfor i = 1:Iteration
        fprintf('index: %d, iter: %d\n', u, i);
        % Initialization
        STAR_Passive_Theta = zeros(M_ite(u),M_ite(u),2);
        MF_Active_Theta = STAR_Passive_Theta;
        STAR_Passive_Theta(:,:,1) = diag(exp(1j*2*pi*rand(M_ite(u),1)));
        STAR_Passive_Theta(:,:,2) = STAR_Passive_Theta(:,:,1);
        MF_Active_Theta(:,:,1) = sqrt(beta_max/2) * STAR_Passive_Theta(:,:,1);
        MF_Active_Theta(:,:,2) = MF_Active_Theta(:,:,1);
        % Position
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Generate_Position(K, User_Position); 
        % Channel
        [h_k,g_k,H] = Generate_Channel(K,M_ite(u),N,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser);
        % passive RIS 
        [R0(u, i)] = Passive_Precoding(N,K,M_ite(u),sigma1,STAR_Passive_Theta,f_k,h_k,g_k,H);
        % active RIS
        [R1(u, i)] = Active_Precoding(N,K,M_ite(u),PRIS,sigma1,sigma2,MF_Active_Theta,f_k,h_k,g_k,H,beta_max);
        % STAR-RIS 
        [R2(u, i)] = STAR_Precoding(N,K,M_ite(u),sigma1,STAR_Passive_Theta,f_k,h_k,g_k,H); 
        % MF-RIS
        [R3(u,i)] = MF_Precoding(N,K,M_ite(u),PRIS,sigma1,sigma2,MF_Active_Theta,f_k,h_k,g_k,H,beta_max);
     end
end
toc;

R0 = mean(R0, 2); 
R1 = mean(R1, 2);
R2 = mean(R2, 2); 
R3 = mean(R3, 2);  
