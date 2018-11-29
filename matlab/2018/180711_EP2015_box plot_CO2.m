%% FC; for INPUT of SDP (seasonal diurnal pattern) box-plot @ SIGMAPLOT
%input = zeros(5142,2);
m = 5142; % length

%% using input (n, 3) << 1 variable + SDP index + (weekday(1) and weekend(0))
n = zeros(192,1); % for temp

output1 = zeros(69,192);
output1 = output1 -999;

for i = 1:m
    n(input(i,2),1) = n(input(i,2),1)+1;
    output1(n(input(i,2),1),input(i,2)) = input(i,1);
end 
clear i

