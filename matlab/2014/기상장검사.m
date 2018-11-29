%%
Tair = nc_varget('Tair_AMIP_2001.nc','Tair');
Tair_temp= zeros(145,192);

for i = 1:2880
    for j = 1:145
        for k = 1:192
            Tair_temp(j,k) = Tair_temp(j,k) + Tair(i,j,k);
        end
    end
end
clear i j k
for j = 1:145
    for k = 1:192
        Tair_temp(j,k) = Tair_temp(j,k)/2880;
    end
end
clear j k

%%
Wind = nc_varget('Wind_AMIP_2001.nc','Wind');
Wind_temp= zeros(145,192);

for i = 1:2880
    for j = 1:145
        for k = 1:192
            Wind_temp(j,k) = Wind_temp(j,k) + Wind(i,j,k);
        end
    end
end
clear i j k
for j = 1:145
    for k = 1:192
        Wind_temp(j,k) = Wind_temp(j,k)/2880;
    end
end
clear j k

%%
SWnet = nc_varget('SWnet_AMIP_2001.nc','SWnet');
SWnet_temp= zeros(145,192);

for i = 1:2880
    for j = 1:145
        for k = 1:192
            SWnet_temp(j,k) = SWnet_temp(j,k) + SWnet(i,j,k);
        end
    end
end
clear i j k
for j = 1:145
    for k = 1:192
        SWnet_temp(j,k) = SWnet_temp(j,k)/2880;
    end
end
clear j k

%%
LWnet = nc_varget('LWnet_AMIP_2001.nc','LWnet');
LWnet_temp= zeros(145,192);

for i = 1:2880
    for j = 1:145
        for k = 1:192
            LWnet_temp(j,k) = LWnet_temp(j,k) + LWnet(i,j,k);
        end
    end
end
clear i j k
for j = 1:145
    for k = 1:192
        LWnet_temp(j,k) = LWnet_temp(j,k)/2880;
    end
end
clear j k

%%
Qair = nc_varget('Qair_AMIP_2001.nc','Qair');
Qair_temp= zeros(145,192);

for i = 1:2880
    for j = 1:145
        for k = 1:192
            Qair_temp(j,k) = Qair_temp(j,k) + Qair(i,j,k);
        end
    end
end
clear i j k
for j = 1:145
    for k = 1:192
        Qair_temp(j,k) = Qair_temp(j,k)/2880;
    end
end
clear j k

%%
PSurf = nc_varget('PSurf_AMIP_2001.nc','PSurf');
PSurf_temp= zeros(145,192);

for i = 1:2880
    for j = 1:145
        for k = 1:192
            PSurf_temp(j,k) = PSurf_temp(j,k) + PSurf(i,j,k);
        end
    end
end
clear i j k
for j = 1:145
    for k = 1:192
        PSurf_temp(j,k) = PSurf_temp(j,k)/2880;
    end
end
clear j k

%%
Rainf = nc_varget('Rainf_AMIP_2001.nc','Rainf');
Rainf_temp= zeros(145,192);

for i = 1:2880
    for j = 1:145
        for k = 1:192
            Rainf_temp(j,k) = Rainf_temp(j,k) + Rainf(i,j,k);
        end
    end
end
clear i j k
for j = 1:145
    for k = 1:192
        Rainf_temp(j,k) = Rainf_temp(j,k)*10800;   %/2880;
    end
end
clear j k