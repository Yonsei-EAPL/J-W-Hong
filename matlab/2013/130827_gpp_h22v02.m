%% extract 
gpp = zeros(240,240,36);
for i = 1:12
    for j = 1:240
        for k = 1:240
            gpp(j,k,i) = gpp2001(j,k,i);
        end
    end
end
for i = 13:24
    for j = 1:240
        for k = 1:240
            gpp(j,k,i) = gpp2002(j,k,i-12);
        end
    end
end
for i = 25:36
    for j = 1:240
        for k = 1:240
            gpp(j,k,i) = gpp2001(j,k,i-24);
        end
    end
end
clear i j k
lai = zeros(240,240,36);
for i = 1:12
    for j = 1:240
        for k = 1:240
            lai(j,k,i) = lai2001(j,k,i);
        end
    end
end
for i = 13:24
    for j = 1:240
        for k = 1:240
            lai(j,k,i) = lai2002(j,k,i-12);
        end
    end
end
for i = 25:36
    for j = 1:240
        for k = 1:240
            lai(j,k,i) = lai2001(j,k,i-24);
        end
    end
end
clear i j k


%% 32767
for i = 1:36
    for j = 1:240
        for k = 1:240
            if gpp(j,k,i)>1000
                gpp(j,k,i) =0;
            end
            if lai(j,k,i)>1000
                lai(j,k,i) =0;
            end
            gpp(j,k,i) = gpp(j,k,i)/10;
            lai(j,k,i) = lai(j,k,i)/10;
        end
    end
end
gpp2 = zeros(36,240,240);
lai2 = zeros(36,240,240);
for i = 1:36
    for j = 1:240
        for k = 1:240
            gpp2(i,j,k) = gpp(k,j,i);
            lai2(i,j,k) = lai(k,j,i);
        end
    end
end
clear i j k


%% 
[lon lat]=tile2geo(22, 2);
lon=double(lon); % convert single to double
lat=double(lat); % convert single to double
lon=imresize(lon, 0.2); % convert 1 km resoultion to 5 km resolution
lat=imresize(lat, 0.2); % convert 1 km resolution to 5 km resoultion
lat2 = zeros(240,1);
lon2 = zeros(240,1);
for i = 1:240
    lat2(i,1) = lat(i,1);
    lon2(i,1) = lon(1,i);
end
clear i lat lon


%% nc file
size_x = 240;
size_y = 240;
size_t = 36;
ncid = netcdf.create('h22v02.nc','NC_WRITE');
lat = netcdf.defDim(ncid,'lat',size_x);
lon = netcdf.defDim(ncid,'lon',size_y);
t = netcdf.defDim(ncid,'t',size_t);
id_gpp = netcdf.defVAR(ncid,'gpp','NC_DOUBLE',[t lon lat]);
id_lai = netcdf.defVAR(ncid,'lai','NC_DOUBLE',[t lon lat]);
id_lat = netcdf.defVAR(ncid,'lat','NC_DOUBLE',[lat]);
id_lon = netcdf.defVAR(ncid,'lon','NC_DOUBLE',[lon]);
netcdf.endDef(ncid);

netcdf.putVar(ncid,id_gpp,gpp2);
netcdf.putVar(ncid,id_lai,lai2);
netcdf.putVar(ncid,id_lat,lat2);
netcdf.putVar(ncid,id_lon,lon2);
netcdf.close(ncid);
clear size_x size_y size_t ncid id_gpp id_lai t lat lon id_lat id_lon



%%
figure(1)
contourf(gpp(:,:,6))

gpp2 = nc_varget('gpp_h22v02.nc','gpp');
%lat2 = nc_varget('gpp_h22v02.nc','lat');
%lon2 = nc_varget('gpp_h22v02.nc','lon');
figure(2)
contourf(gpp2(:,:,6))
