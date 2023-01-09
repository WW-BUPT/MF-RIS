function [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Generate_Position(K, user_X)
Dis_BStoUser=zeros(K,1);
Dis_RIStoUser=zeros(K,1);

BS_position = [0 0 0];
RIS_position = [30 50 0];

user_position = zeros(K,3);
user_position(1,:)=[user_X 40 0];
user_position(2,:)=[user_X 60 0]; 

Dis_BStoRIS = Generate_Distance(BS_position,RIS_position);

for k=1:K
    user_position_temp = reshape(user_position(k,:),3,1);
    Dis_BStoUser(k) = Generate_Distance(BS_position,user_position_temp);
end

for k=1:K
    user_position_temp = reshape(user_position(k,:),3,1);
    Dis_RIStoUser(k) = Generate_Distance(RIS_position,user_position_temp);
end

end