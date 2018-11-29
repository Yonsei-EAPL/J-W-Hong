ah = nc_varget('AH.nc','AH');
lat = nc_varget('AH.nc','LAT');
lon = nc_varget('AH.nc','LON');

%EP  37.635033¡ÆN 126.928661¡ÆE

ep_lat = 37.635033;
ep_lon = 126.928661;

for i = 1:205
    for j = 1:215
        lat(i,j) = lat(i,j)-ep_lat;
        lon(i,j) = lon(i,j)-ep_lon;
    end
end
clear i j

dis = zeros(205,215);
for i = 1:205
    for j = 1:215
        dis(i,j) = lat(i,j)^2 + lon(i,j)^2;
    end
end
clear i j

lat = nc_varget('AH.nc','LAT');
lon = nc_varget('AH.nc','LON');

ah108 = ah(1:24,111,130);
%ahep = ah(1:24,x,y);
ahep = ah(1:24,117,126);


% for i = 1:10
%     for j = 1:10
%         k = zeros(24,1);
%         k = ah(1:24,x-5+i,y-5+j);
%         hold on 
%         plot(k,'-r')
%     end
% end
% clear i j
% 
ah24 = zeros(205,215);
for i = 1:205
    for j = 1:215
        for k = 1:24
            ah24(i,j) = ah24(i,j)+ah(k,i,j); 
        end
        ah24(i,j) = ah24(i,j)/24;
    end
end
clear i j k
        

x=zeros(25,1);
y=x;
for i = 1:25
    x(i,1) = (i-1)*2;
    if i==1
        y(i,1) = ahep(24,1);
    else
        y(i,1) = ahep(i-1,1);
    end
end 
clear i

x2 = zeros(49,1);
for i = 1:49
    x2(i,1) = i-1;
end
clear i

y2 = interp1(x,y,x2);





