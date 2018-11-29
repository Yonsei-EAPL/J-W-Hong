% data_108 (hourly)
% 1: YYYYMMDDHH
% 2: YYYYMMDD00
% 3: rainfall
% 4: solar radiation (MJ m-2)
% 5: rainfall (1/0)
% 6: hour after rainfall
% 7: solar radiation after rainfall

% data_108_2 (daily)
% 1: YYYYMMDDHH
% 2: rainfall (mm day-1)
% 3: rainfall (1/0)

% data_UHI
% 1: YYYYMMDD00
% 2: site number
% 3: daily maximum UHI (T_AWS - T_±èÆ÷°øÇ×)
% 4: LCZ at maximum UHI
% 5: HH (time)
% 6: number of all data (24/23 are OK)
% 7: number of good data
% 8: number of gaps
% 9: rainy day(1/0)
% 10: hour after rainfall
% 11: solar radiation after rainfall

data_uhi(1,9)=0;
data_uhi(1,10)=0;
data_uhi(1,11)=0;
for i = 1:length(data_uhi(:,1))
    a = find(data_108_2(:,1)==data_uhi(i,1));
    data_uhi(i,9) = data_108_2(a,3);
    if data_uhi(i,5)>20000
        a = find(data_108(:,1)==data_uhi(i,5));
        data_uhi(i,10) = data_108(a,6);
        data_uhi(i,11) = data_108(a,7);    
        clear a    
    end
end
clear i

save('data_uhi.mat','data_uhi');
