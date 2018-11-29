% data Ã³¸®1
% 1: site number
% 2: winddirection
% 3: LCZ
% 4: YYYYMMDDHH
% 5: air temperature
% 6: YYYYMMDD00 (³¯Â¥ - for ÇÇ¹ş)
% 7: T_AWS - T_±èÆ÷°øÇ×

load('data_hour.mat');
load('data_kp.mat');

data_hour(1,7) = 0;
for i = 1:length(data_hour)
    if (data_hour(i,5)>-40)&&(data_hour(i,5)<45)
        temp = find(data_kp(:,1)==data_hour(i,4));
        if (data_kp(temp,2)>-40)&&(data_kp(temp,2)<45)
            data_hour(i,7) = data_hour(i,5)- data_kp(temp,2);
        else
            data_hour(i,7) = 888.8;
        end
    end    
end
clear i temp

save('data_hour2.mat','data_hour');



