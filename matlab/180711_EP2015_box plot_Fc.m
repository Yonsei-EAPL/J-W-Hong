%% FC; for INPUT of SDP (seasonal diurnal pattern) box-plot @ SIGMAPLOT
%input = zeros(5002,3);
m = 5002; % length

%% using input (n, 3) << 1 variable + SDP index + (weekday(1) and weekend(0))
n = zeros(192,3); % for temp ; total; weekday; weekend

output1 = zeros(69,192);
output2 = zeros(49,192);
output3 = zeros(22,192);
output1 = output1 -999;
output2 = output2 -999;
output3 = output3 -999;

for i = 1:m
    if input(i,3)==1 %weekday
        n(input(i,2),1) = n(input(i,2),1)+1;
        n(input(i,2),2) = n(input(i,2),2)+1;
        output1(n(input(i,2),1),input(i,2)) = input(i,1);
        output2(n(input(i,2),2),input(i,2)) = input(i,1);
    else %weekend
        n(input(i,2),1) = n(input(i,2),1)+1;
        n(input(i,2),3) = n(input(i,2),3)+1;
        output1(n(input(i,2),1),input(i,2)) = input(i,1);
        output3(n(input(i,2),3),input(i,2)) = input(i,1);
    end
end 
clear i

