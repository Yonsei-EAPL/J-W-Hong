% data 처리1
% 1: site number
% 2: winddirection
% 3: LCZ
% 4: YYYYMMDDHH
% 5: air temperature
% 6: YYYYMMDD00 (날짜 - for 피벗)

load('data_hour.mat');

data_hour(1,6) = 0;
for i = 1:length(data_hour)
    data_hour(i,6) = data_hour(i,4)- mod(data_hour(i,4),100);
end
clear i

save('data_hour.mat','data_hour');



