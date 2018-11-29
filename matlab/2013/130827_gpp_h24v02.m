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


%% nc file

size_x = 240;
size_y = 240;
size_t = 36;

ncid = netcdf.create('gpp_h24v02.nc','NC_WRITE');

lat = netcdf.defDim(ncid,'lat',size_x);
lon = netcdf.defDim(ncid,'lon',size_y);
t = netcdf.defDim(ncid,'time',size_t);

gpp_h24v02 = netcdf.defVAR(ncid,'gpp','NC_FLOAT',[t lon lat]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,gpp_h24v02,gpp);
netcdf.close(ncid);

clear size_x size_y size_t ncid lon lat gpp_h24v02 t