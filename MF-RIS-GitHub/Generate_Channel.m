function [h_k,g_k,H] = Generate_Channel(K,M,N,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser)
h_k=zeros(K,N);
g_k=zeros(K,M);
H=zeros(M,N);

for k=1:K
	h_k(k,:)=Generate_Channel_h(N,1,Dis_BStoUser(k));             
end


for k=1:K
	g_k(k,:)=Generate_Channel_g(M,1,Dis_RIStoUser(k));
end

H(:,:)=Generate_Channel_HH(M,N,Dis_BStoRIS);

end

