% data_108 (hourly)
% 1: YYYYMMDDHH
% 2: YYYYMMDD00
% 3: rainfall
% 4: solar radiation (MJ m-2)
% 5: rainfall (1/0)
% 6: hour after rainfall
% 7: solar radiation after rainfall

for i =1:length(data_108(:,1))
    if data_108(i,3)>0
        data_108(i,5) = 1;
        data_108(i,6) = 0;
        data_108(i,7) = 0;
    else
        if i >1
            data_108(i,5) = 0;
            data_108(i,6) = data_108(i-1,6)+1;
            data_108(i,7) = data_108(i-1,7)+data_108(i,4);
        end
    end
end
clear i

% data_108_2 (daily)
% 1: YYYYMMDDHH
% 2: rainfall (mm day-1)
% 3: rainfall (1/0)
data_108_2 = unique(data_108(:,2));
data_108_2(1,2) =0;
for i = 1:length(data_108(:,1))
    if data_108(i,3)>0
        a = find(data_108_2(:,1)==data_108(i,2));
        data_108_2(a,2) = data_108_2(a,2)+data_108(i,3);
    end
    clear a
end
clear i
data_108_2(1,3) =0;
for i = 1:length(data_108_2(:,1))
    if data_108_2(i,2)>0
        data_108_2(i,3) = 1;
    end
end
clear i

save('data_108.mat','x','data_108','data_108_2');