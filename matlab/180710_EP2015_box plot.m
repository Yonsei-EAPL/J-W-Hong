%% E-flux; for INPUT of SDP (seasonal diurnal pattern) box-plot @ SIGMAPLOT
% input = zeros(4618,9);
m  = 4618; % length

%% using input (n, 9) << 8 variables + SDP index 
n = zeros(192,1); % for temp

output = zeros(8,66,192); % total output
output1 = zeros(66,192); %QH
output2 = zeros(66,192); %QE 
output3 = zeros(66,192); %Kdn
output4 = zeros(66,192); %Kup
output5 = zeros(66,192); %Ldn
output6 = zeros(66,192); %Lup
output7 = zeros(66,192); %Q*
output8 = zeros(66,192); %residual

output = output -999;
output1 = output1 -999;
output2 = output2 -999;
output3 = output3 -999;
output4 = output4 -999;
output5 = output5 -999;
output6 = output6 -999;
output7 = output7 -999;
output8 = output8 -999;

for i = 1:m
    n(input(i,9),1) = n(input(i,9),1)+1;
    for j = 1:8
        output(j,n(input(i,9),1),input(i,9))=input(i,j);
    end
    clear j
end
clear i

for i= 1:66
    for j = 1:192
        output1(i,j) = output(1,i,j);
        output2(i,j) = output(2,i,j);
        output3(i,j) = output(3,i,j);
        output4(i,j) = output(4,i,j);
        output5(i,j) = output(5,i,j);
        output6(i,j) = output(6,i,j);
        output7(i,j) = output(7,i,j);
        output8(i,j) = output(8,i,j);
    end
end
clear i j 
    


