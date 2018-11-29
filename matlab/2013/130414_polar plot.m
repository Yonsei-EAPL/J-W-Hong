temp_max = 8;
temp_min = 0;


WD_backup = WD;
data_backup = data;
a = size(WD_backup);
sizeb = a(1,1);
clear a
for i = 1:sizeb
    if WD_backup(i,1)*(-1)+90<0
        WD_backup(i,1) = WD_backup(i,1)*(-1)+450;
    else
        WD_backup(i,1) = WD_backup(i,1)*(-1)+90;
    end
end
% for i = 1:sizeb
%     data_backup(i,1) = log10(data_backup(i,1));
% end
for i = 1:sizeb
    if data_backup(i,1) > temp_max
        data_backup(i,1) =temp_max;
    elseif data_backup(i,1) < temp_min
        data_backup(i,1) = temp_min;
    end
end

WD_backup2 = WD2;
data_backup2 = data2;
a = size(WD_backup2);
sizeb = a(1,1);
clear a
for i = 1:sizeb
    if WD_backup2(i,1)*(-1)+90<0
        WD_backup2(i,1) = WD_backup2(i,1)*(-1)+450;
    else
        WD_backup2(i,1) = WD_backup2(i,1)*(-1)+90;
    end
end
for i = 1:sizeb
    if data_backup2(i,1) > temp_max
        data_backup2(i,1) =temp_max;
    elseif data_backup2(i,1) < temp_min
        data_backup2(i,1) =temp_min;
    end
end
% for i = 1:sizeb
%     data_backup2(i,1) = log10(data_backup2(i,1));
% end
% 
% max_1 = max(data_backup(:,1));
% max_2 = max(data_backup2(:,1));
% 
% 
% figure(1)
% wind_rose(WD_backup(:,1),data_backup(:,1))
% figure(2)
% wind_rose(WD_backup2(:,1),data_backup2(:,1))
[max_1, max_1_p] = max(data_backup(:,1));
[max_2, max_2_p] = max(data_backup2(:,1));
if max_1>max_2
    data_backup2(max_2_p,1) = max_1;
else
    data_backup(max_1_p,1) = max_2;
end
[max_1, max_1_p] = min(data_backup(:,1));
[max_2, max_2_p] = min(data_backup2(:,1));
if max_1>max_2
    data_backup(max_1_p,1) = max_2;
else
    data_backup2(max_2_p,1) = max_1;
end
figure(1)
% wind_rose(WD_backup(:,1),data_backup(:,1))
wind_rose(WD_backup(:,1),data_backup(:,1),'ci',[5,10,15,20,25])
% wind_rose(WD_backup(:,1),data_backup(:,1),'ci',[2,4,6,8])
% wind_rose(WD_backup(:,1),data_backup(:,1),'ci',[3,6,9,12,15])
figure(2)
% wind_rose(WD_backup2(:,1),data_backup2(:,1))
% wind_rose(WD_backup2(:,1),data_backup2(:,1),'ci',[3,6,9,12,15])
% wind_rose(WD_backup2(:,1),data_backup2(:,1),'ci',[2,4,6,8])
wind_rose(WD_backup2(:,1),data_backup2(:,1),'ci',[5,10,15,20,25])
clear i max_1 max_2 max_1_p max_2_p sizeb 



% for i = 1:sizeb
%     a = WD(i,1) + 5;
% 	if a>360
%         WD(i,1) = WD(i,1) + 5 - 360;
%     end
%     
% end
% 
% temp = zeros(36,2);
% for i = 1:sizeb
%     WDi = fix(WD(i,1)/10+1);
%     temp(WDi,1) = temp(WDi,1)+1;
%     temp(WDi,2) = temp(WDi,2)+data(i,1);
% end
% for i = 1:36
%     temp(i,2) = temp(i,2)/temp(i,1);
% end
% for i = 1:36
%     temp(i,1) = (i*10-5)/180*pi();
% end
% clear i 




