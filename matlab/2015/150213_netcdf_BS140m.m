[a b]=size(data);

for i = 1:a
    Ui = data(:,1);
    WDi = data(:,2);
    ustari = data(:,3);
    Tsi = data(:,4);
    covwTsi = data(:,5);
end
clear i

size_time = a;
ncid = netcdf.create('BS140m.nc','NC_WRITE');

time = netcdf.defDim(ncid,'time',size_time);

idtime = netcdf.defVAR(ncid,'time','NC_FLOAT',[time]);
netcdf.putAtt(ncid,idtime,'units','yyyyMMddhhmm');
idU = netcdf.defVAR(ncid,'U','NC_FLOAT',[time]);
netcdf.putAtt(ncid,idU,'units','m');
idWD = netcdf.defVAR(ncid,'WD','NC_FLOAT',[time]);
netcdf.putAtt(ncid,idWD,'units','degree');
idustar = netcdf.defVAR(ncid,'ustar','NC_FLOAT',[time]);
netcdf.putAtt(ncid,idustar,'units','m s-1');
idTs = netcdf.defVAR(ncid,'Ts','NC_FLOAT',[time]);
netcdf.putAtt(ncid,idTs,'units','degreeC');
idcovwTs = netcdf.defVAR(ncid,'covwTs','NC_FLOAT',[time]);
netcdf.putAtt(ncid,idcovwTs,'units','degreeC m s-1');
netcdf.endDef(ncid);

netcdf.putVar(ncid,idtime,timei);
netcdf.putVar(ncid,idU,Ui);
netcdf.putVar(ncid,idWD,WDi);
netcdf.putVar(ncid,idustar,ustari);
netcdf.putVar(ncid,idTs,Tsi);
netcdf.putVar(ncid,idcovwTs,covwTsi);
netcdf.close(ncid);