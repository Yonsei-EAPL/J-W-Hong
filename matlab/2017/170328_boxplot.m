%% input description
% data
% 1 : season (from 1 (spring, MAM) to 4 (winter, DJF))
% 2 : day_48 (30min interval, from 1 (00:15) to 48 (23:45))
% 3-10 : QH, QE, Kdn, Kup, Ldn, Lup, Q*, dQS

%% output description
% output1 (QH) : maximum number * 48 (day48) * 4 (season) 
% 2-8 : QE, Kdn, Kup, Ldn, Lup, Q*, dQS

%% statistics
[n_data n_var] = size(data);

n_statistics = zeros(48,4);
for i = 1:n_data
    n_statistics(data(i,2),data(i,1)) = n_statistics(data(i,2),data(i,1))+1;
end
clear i

output1 = zeros(max(max(n_statistics(:,:))),48*4);
output2 = zeros(max(max(n_statistics(:,:))),48*4);
output3 = zeros(max(max(n_statistics(:,:))),48*4);
output4 = zeros(max(max(n_statistics(:,:))),48*4);
output5 = zeros(max(max(n_statistics(:,:))),48*4);
output6 = zeros(max(max(n_statistics(:,:))),48*4);
output7 = zeros(max(max(n_statistics(:,:))),48*4);
output8 = zeros(max(max(n_statistics(:,:))),48*4);

temp = zeros(48,4);
for i = 1:n_data
    temp(data(i,2),data(i,1)) = temp(data(i,2),data(i,1))+1;
    output1(temp(data(i,2),data(i,1)),(data(i,1)-1)*48+data(i,2)) = data(i,3);
    output2(temp(data(i,2),data(i,1)),(data(i,1)-1)*48+data(i,2)) = data(i,4);
    output3(temp(data(i,2),data(i,1)),(data(i,1)-1)*48+data(i,2)) = data(i,5);
    output4(temp(data(i,2),data(i,1)),(data(i,1)-1)*48+data(i,2)) = data(i,6);
    output5(temp(data(i,2),data(i,1)),(data(i,1)-1)*48+data(i,2)) = data(i,7);
    output6(temp(data(i,2),data(i,1)),(data(i,1)-1)*48+data(i,2)) = data(i,8);
    output7(temp(data(i,2),data(i,1)),(data(i,1)-1)*48+data(i,2)) = data(i,9);
    output8(temp(data(i,2),data(i,1)),(data(i,1)-1)*48+data(i,2)) = data(i,10);
end
clear i

for i = 1:max(n_statistics(:,:))
    for j = 1:48*4
        if output1(i,j)==0
            output1(i,j)=-999;
            output2(i,j)=-999;
            output3(i,j)=-999;
            output4(i,j)=-999;
            output5(i,j)=-999;
            output6(i,j)=-999;
            output7(i,j)=-999;
            output8(i,j)=-999;
        end
    end
end
clear i j