timei = nc_varget('BSground.nc','time');    

[a b]=size(data);

for i = 1:a
	if data(i,4)~=-999
		data(i,4) = data(i,4)-273.15;
	end
end
clear i

% for i = 1:a
    Ui = data(:,1);
    WDi = data(:,2);
    ustari = data(:,3);
    Tsi = data(:,4);
    covwTsi = data(:,5);
% end
% clear i

size_time = a;
ncid = netcdf.create('BS140m.nc','NC_WRITE');

time = netcdf.defDim(ncid,'time',size_time);

idtime = netcdf.defVar(ncid,'time','NC_DOUBLE',[time]);
netcdf.putAtt(ncid,idtime,'units','yyyyMMddhhmm');
idU = netcdf.defVar(ncid,'U','NC_DOUBLE',[time]);
netcdf.putAtt(ncid,idU,'units','m s-1');
idWD = netcdf.defVar(ncid,'WD','NC_DOUBLE',[time]);
netcdf.putAtt(ncid,idWD,'units','degree');
idustar = netcdf.defVar(ncid,'ustar','NC_DOUBLE',[time]);
netcdf.putAtt(ncid,idustar,'units','m s-1');
idTs = netcdf.defVar(ncid,'Ts','NC_DOUBLE',[time]);
netcdf.putAtt(ncid,idTs,'units','degreeC');
idcovwTs = netcdf.defVar(ncid,'covwTs','NC_DOUBLE',[time]);
netcdf.putAtt(ncid,idcovwTs,'units','degreeC m s-1');
netcdf.endDef(ncid);

netcdf.putVar(ncid,idtime,timei);
netcdf.putVar(ncid,idU,Ui);
netcdf.putVar(ncid,idWD,WDi);
netcdf.putVar(ncid,idustar,ustari);
netcdf.putVar(ncid,idTs,Tsi);
netcdf.putVar(ncid,idcovwTs,covwTsi);
netcdf.close(ncid);