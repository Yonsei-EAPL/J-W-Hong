sat_up = zeros(1800,1800,7);
for m = 1:7
    for i = 1:1800
        for j = 1:1800
            for k = 1:2
                for l = 1:2
                    sat_up(i,j,m) = sat_up(i,j,m)+ sat_urb((i-1)*2+k,(j-1)*2+l,m);
                end
            end
            sat_up(i,j,m) = sat_up(i,j,m)/4;
        end
    end
end
clear i j k l m

dmsp_urb = zeros(1800,1800,31);
for i = 1:1800
    for j = 1:1800
        dmsp_urb(i,j,1) = f101992_land(i,j);
        dmsp_urb(i,j,2) = f101993_land(i,j);
        dmsp_urb(i,j,3) = f101994_land(i,j);
        dmsp_urb(i,j,4) = f121994_land(i,j);
        dmsp_urb(i,j,5) = f121995_land(i,j);
        dmsp_urb(i,j,6) = f121996_land(i,j);
        dmsp_urb(i,j,7) = f121997_land(i,j);
        dmsp_urb(i,j,8) = f121998_land(i,j);
        dmsp_urb(i,j,9) = f121999_land(i,j);
        dmsp_urb(i,j,10) = f141997_land(i,j);
        dmsp_urb(i,j,11) = f141998_land(i,j);
        dmsp_urb(i,j,12) = f141999_land(i,j);
        dmsp_urb(i,j,13) = f142000_land(i,j);
        dmsp_urb(i,j,14) = f142001_land(i,j);
        dmsp_urb(i,j,15) = f142002_land(i,j);
        dmsp_urb(i,j,16) = f142003_land(i,j);
        dmsp_urb(i,j,17) = f152000_land(i,j);
        dmsp_urb(i,j,18) = f152001_land(i,j);
        dmsp_urb(i,j,19) = f152002_land(i,j);
        dmsp_urb(i,j,20) = f152003_land(i,j);
        dmsp_urb(i,j,21) = f152004_land(i,j);
        dmsp_urb(i,j,22) = f152005_land(i,j);
        dmsp_urb(i,j,23) = f152006_land(i,j);
        dmsp_urb(i,j,24) = f152007_land(i,j);
        dmsp_urb(i,j,25) = f162004_land(i,j);
        dmsp_urb(i,j,26) = f162005_land(i,j);
        dmsp_urb(i,j,27) = f162006_land(i,j);
        dmsp_urb(i,j,28) = f162007_land(i,j);
        dmsp_urb(i,j,29) = f162008_land(i,j);
        dmsp_urb(i,j,30) = f162009_land(i,j);
        dmsp_urb(i,j,31) = f182010_land(i,j);
    end
end
clear i j

dmsp_seoul = zeros(150,150,31);
for i = 1:31
    for j = 1:150
        for k = 1:150
            dmsp_seoul(j,k,i) = dmsp_urb(850+j,750+k,i);
        end
    end
end
clear i j k 

sat_seoul = zeros(150,150,7);
for i = 1:7
    for j = 1:150
        for k = 1:150
            sat_seoul(j,k,i) = sat_up(850+j,750+k,i);
        end
    end
end
clear i j k 

data_mask_seoul = zeros(150,150);
for i = 1:150
    for j = 1:150
        data_mask_seoul(i,j) = data_mask(850+i,750+j);
    end
end
clear i j 


%% calculation the threshold for seoul_raw(1-63)
seoul_th_dmsp = zeros(63,31);
for j = 1:31
    temp = zeros(150,150);
    for k = 1:150
        for l = 1:150
            temp(k,l) = seoul(k,l,j);
        end
    end
    for i = 1:63
        for k = 1:150
            for l = 1:150
                if temp(k,l)<i
                    temp(k,l) = 0;
                end
            end
        end
        [B, L] = bwboundaries(temp,'holes');
        for k = 1:length(B)
            boundary = B{k};
            delta_sq = diff(boundary).^2;
            perimeter = sum(sqrt(sum(delta_sq,2)));
            seoul_th_dmsp(i,j) = seoul_th_dmsp(i,j)+perimeter;
        end
    end
end
clear i j k l B L boundary delta_sq perimeter 

th = 61.5;
dmsp_urb = dmsp_urb_raw;
for i = 1:31
    for j = 1:1800
        for k = 1:1800
            if dmsp_urb(j,k,i)<th
                dmsp_urb(j,k,i)=0;
            end                
        end
    end
end
clear i j k

%% filter satellite data (delete 25%)
for i = 1:7
    for j = 1:1800
        for k = 1:1800
            if sat_up(j,k,i)==25
                sat_up(j,k,i) = 0;
            end
        end
    end
end
clear i j k


%% thresold from satellite (using all)
n_sat = zeros(7,1);
for i = 1:7
    temp = 0;
    for j = 1:1800
        for k = 1:1800
            if sat_up(j,k,i)>0
                temp = temp+1;
            end
        end
    end
    n_sat(i,1) = temp;
end
clear i j k temp

th_sat = zeros(63,14);
num_dmsp = [14,18;15,19;16,20;21,25;22,26;23,27;24,28];
for i = 1:7
    for j = 1:2
        for k = 1:63
            x1=0;
            x2=0;
            for l = 1:1800
                for m = 1:1800
%                     if (sat_up(l,m,i)>0)&&(dmsp_urb(l,m,num_dmsp(i,j))>=k)
%                         x1 = x1+1;
%                     end
%                     if (sat_up(l,m,i)==0)&&(dmsp_urb(l,m,num_dmsp(i,j))>=k)
%                         x2 = x2+1;
%                     end
                    if dmsp_urb(l,m,num_dmsp(i,j))>=k
                        x1 = x1 + 1;
                    end
                end
            end
            th_sat(k,(i-1)*2+j)= (x1)/n_sat(i,1);
        end
    end
end
clear i j k l m x1 x2 


%% thresold from satellite (using seoul)
n_sat = zeros(7,1);
for i = 1:7
    temp = 0;
    for j = 1:150
        for k = 1:150
            if sat_seoul(j,k,i)>0
                temp = temp+1;
            end
        end
    end
    n_sat(i,1) = temp;
end
clear i j k temp

th_sat = zeros(63,14);
num_dmsp = [14,18;15,19;16,20;21,25;22,26;23,27;24,28];
for i = 1:7
    for j = 1:2
        for k = 1:63
            x1=0;
            x2=0;
            for l = 1:150
                for m = 1:150
%                     if (sat_seoul(l,m,i)>0)&&(dmsp_seoul(l,m,num_dmsp(i,j))>=k)
%                         x1 = x1+1;
%                     end
%                     if (sat_seoul(l,m,i)==0)&&(dmsp_seoul(l,m,num_dmsp(i,j))>=k)
%                         x2 = x2+1;
%                     end
                    if dmsp_seoul(l,m,num_dmsp(i,j))>=k
                        x1 = x1+1;
                    end
                end
            end
            th_sat(k,(i-1)*2+j)= (x1)/n_sat(i,1);
        end
    end
end
clear i j k l m x1 x2 

for i = 1:63
    for j = 1:14
        th_sat(i,j) = abs(1-th_sat(i,j));
    end
end
clear i j

th_sat_final = zeros(14,1);
for i = 1:14
    [a b] = min(th_sat(:,i));
    th_sat_final(i,1) = b;
end
clear i


%% averaging two year 2001-2007
avg_sat = zeros(1800,1800,7);
for i = 1:7
    for j = 1:2
        for k = 1:1800
            for l = 1:1800
                avg_sat(k,l,i) = avg_sat(k,l,i)+dmsp_urb(k,l,num_dmsp(i,j))/2;
            end
        end
    end
end
clear i j k l

th_avg_sat = zeros(63,7);
for i = 1:7
    for j = 1:63
        x1 = 0;
        for k = 1:1800
            for l = 1:1800
            end
        end
    end
end
clear i j k l






%% 그림 그릴 때,
year = 7;
temp = zeros(150,150);
for i = 1:150
    for j = 1:150
        temp(i,j) = seoul(i,j,year);
    end
end
clear i j year
image(temp);
hold on
spy(data_mask_seoul,'k');

year = 19;
temp = zeros(1800,1800);
for i = 1:1800
    for j = 1:1800
        if dmsp_urb(i,j,year)>=40
            temp(i,j) = dmsp_urb(i,j,year);
        end
    end
end
clear i j year
image(temp);
hold on
spy(data_mask,'k');

sat = 2;
temp = zeros(1800,1800);
for i = 1:1800
    for j = 1:1800
        temp(i,j) = sat_up(i,j,sat);
    end
end
clear i j sat
image(temp);
hold on
spy(data_mask,'k');

