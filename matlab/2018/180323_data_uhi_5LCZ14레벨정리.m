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
% 12: UHI final with environmental lapse rate
% 13: LCZ-14levels

data_uhi(1,13)=0;
for i = 1:length(data_uhi(:,1))
    if (data_uhi(i,4)>0)&&(data_uhi(i,4)<300)
        data_uhi(i,13) = list_lcz(data_uhi(i,4),2);
    end
end
clear i

save('data_uhi.mat','data_uhi');
