Qe_u0 = nc_varget('gdk.tstep.nc','latent_heat');
Qh_u0 = nc_varget('gdk.tstep.nc','ftl_gb');

Qe_u50 = nc_varget('gdk.tstep.nc','latent_heat');
Qh_u50 = nc_varget('gdk.tstep.nc','ftl_gb');

Qe_u90 = nc_varget('gdk.tstep.nc','latent_heat');
Qh_u90 = nc_varget('gdk.tstep.nc','ftl_gb');

Qe_u100 = nc_varget('gdk.tstep.nc','latent_heat');
Qh_u100 = nc_varget('gdk.tstep.nc','ftl_gb');
Qa = nc_varget('gdk.tstep.nc','anthrop_heat');



data = zeros(52608,9);
for i = 1:52608
    data(i,1) = Qa(i,6);
    data(i,2) = Qh_u100(i,1);
    data(i,3) = Qe_u100(i,1);
    data(i,4) = Qh_u90(i,1);
    data(i,5) = Qe_u90(i,1);
    data(i,6) = Qh_u50(i,1);
    data(i,7) = Qe_u50(i,1);
    data(i,8) = Qh_u0(i,1);
    data(i,9) = Qe_u0(i,1);
end
clear i