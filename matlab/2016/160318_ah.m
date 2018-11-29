ah = nc_varget('AH.nc','AH');
ah_lat = nc_varget('AH.nc','LAT');
ah_lon = nc_varget('AH.nc','LON');

ah_mean = zeros(205,215);
for i = 1:205
    for j = 1:215
        ah_mean(i,j) = mean(ah(:,i,j));
    end
end
clear i j 

ah_mean_list = zeros(205*215,1);
for i = 1:205
    for j = 1:215
        ah_mean_list((i-1)*215+j,1) = ah_mean(i,j);
    end
end
clear i j
[x y] = find(ah_mean_list(:,1)>0);
ah_mean_list_azero = zeros(length(y),1);
for i = 1:length(y)
    ah_mean_list_azero(i,1) = ah_mean_list(x(i,1),1);
end
clear i
hist(ah_mean_list_azero,0:2:65)

FRC2 = zeros(205,215);
for i =1:215
    for j = 1:205
        FRC2(j,i) = FRC(i,j);
    end
end
clear i j

FRC3 = FRC2;
ah_mean2 = ah_mean;
for i = 1:205
    for j = 1:215
        if (ah_mean(i,j)==0)&&(FRC2(i,j)>0)
            FRC3(i,j)=0;
        elseif (ah_mean(i,j)>0)&&(FRC2(i,j)==0)
            ah_mean(i,j)=0;
        end
    end
end
clear i j

ah_mean_list = zeros(205*215,1);
frc_list = zeros(205*215,1);
for i = 1:205
    for j = 1:215
        ah_mean_list((i-1)*215+j,1) = ah_mean(i,j);
        frc_list((i-1)*215+j,1) = FRC3(i,j);
    end
end
clear i j
[x y] = find(ah_mean_list(:,1)>0);
ah_mean_list_azero = zeros(length(y),1);
frc_list_azero = zeros(length(y),1);
for i = 1:length(y)
    ah_mean_list_azero(i,1) = ah_mean_list(x(i,1),1);
    frc_list_azero(i,1) = frc_list(x(i,1),1);
end
clear i




