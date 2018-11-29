%% process the BLH data from ECMWF reanalysis
% Final request
% 
% Stream:	Atmospheric model
% Area:	37.7¡ÆN 126.77¡ÆE 37.4¡ÆN 127.16¡ÆE
% Parameter:	Boundary layer height
% Dataset:	interim_daily
% Step:	3 to 12 by 3
% Version:	1
% Type of level:	Surface
% Time:	00:00:00, 12:00:00
% Date:	20010101 to 20011231, 20020201 to 20160229
% Grid:	0.25¡Æ x 0.25¡Æ
% Type:	Forecast
% Class:	ERA Interim

blh_seoul = zeros(44056,1);
blh_seoul(:,1) = blh(:,2,2);

time = nc_varget('_grib2netcdf-atls15-95e2cf679cd58ee9b4db4dd119a05a8d-Ok0tay.nc','time');
blh = nc_varget('_grib2netcdf-atls15-95e2cf679cd58ee9b4db4dd119a05a8d-Ok0tay.nc','blh');
time2 = double(time)/24;

