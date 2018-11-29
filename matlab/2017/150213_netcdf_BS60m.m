timei = nc_varget('BSground.nc','time');    

[a b]=size(data);
% 
%     
%     for i = 1:a
%     data(i,15) = data(i,11)-data(i,12)+data(i,13)-data(i,14);
%     data(i,4) = data(i,4)-273.15;
% end
% clear i

for i = 1:a
	if data(i,4)~=-999
		data(i,4) = data(i,4)-273.15;
	end
end
clear i

%     for i = 1:a
        Ui = data(:,1);
        WDi = data(:,2);
        ustari = data(:,3);
        Tsi = data(:,4);
        CO2i = data(:,5);
        H2Oi = data(:,6);
        covwTsi = data(:,7);
        Qhi = data(:,8);
        Qei = data(:,9);
        Fci = data(:,10);
        Kdni = data(:,11);
        Kupi = data(:,12);
        Ldni = data(:,13);
        Lupi = data(:,14);
        Qstari = data(:,15);
%     end
%     clear i

    size_time = a;
    ncid = netcdf.create('BS60m.nc','NC_WRITE');

    time = netcdf.defDim(ncid,'time',size_time);
    
    idtime = netcdf.defVar(ncid,'time','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idtime,'units','yyyyMMddhhmm');
    idU = netcdf.defVar(ncid,'U','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idU,'units','m');
    idWD = netcdf.defVar(ncid,'WD','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idWD,'units','degree');
    idustar = netcdf.defVar(ncid,'ustar','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idustar,'units','m s-1');
    idTs = netcdf.defVar(ncid,'Ts','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idTs,'units','degreeC');
    idCO2 = netcdf.defVar(ncid,'CO2','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idCO2,'units','ppm');
    idH2O = netcdf.defVar(ncid,'H2O','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idH2O,'units','ppt');
    idcovwTs = netcdf.defVar(ncid,'covwTs','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idcovwTs,'units','degreeC m s-1');
    idQh = netcdf.defVar(ncid,'Qh','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idQh,'units','W m-2');
    idQe = netcdf.defVar(ncid,'Qe','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idQe,'units','W m-2');
    idFc = netcdf.defVar(ncid,'Fc','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idFc,'units','umol m-2 s-1');
    idKdn = netcdf.defVar(ncid,'Kdn','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idKdn,'units','W m-2');
    idKup = netcdf.defVar(ncid,'Kup','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idKup,'units','W m-2');
    idLdn = netcdf.defVar(ncid,'Ldn','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idLdn,'units','W m-2');
    idLup = netcdf.defVar(ncid,'Lup','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idLup,'units','W m-2');
    idQstar = netcdf.defVar(ncid,'Qstar','NC_DOUBLE',[time]);
    netcdf.putAtt(ncid,idQstar,'units','W m-2');
    netcdf.endDef(ncid);

    netcdf.putVar(ncid,idtime,timei);
    netcdf.putVar(ncid,idU,Ui);
    netcdf.putVar(ncid,idWD,WDi);
    netcdf.putVar(ncid,idustar,ustari);
    netcdf.putVar(ncid,idTs,Tsi);
    netcdf.putVar(ncid,idCO2,CO2i);
    netcdf.putVar(ncid,idH2O,H2Oi);
    netcdf.putVar(ncid,idcovwTs,covwTsi);
    netcdf.putVar(ncid,idQh,Qhi);
    netcdf.putVar(ncid,idQe,Qei);
    netcdf.putVar(ncid,idFc,Fci);
    netcdf.putVar(ncid,idKdn,Kdni);
    netcdf.putVar(ncid,idKup,Kupi);
    netcdf.putVar(ncid,idLdn,Ldni);
    netcdf.putVar(ncid,idLup,Lupi);
    netcdf.putVar(ncid,idQstar,Qstari);
    netcdf.close(ncid);