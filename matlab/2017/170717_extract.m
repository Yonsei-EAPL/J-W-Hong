t_j = nc_varget('jc.nc','time'); % 1407231530 부터 30분 간격
t_m = nc_varget('mc.nc','Times'); % 1407231500 부터 15분 간격

rn_jc = nc_varget('jc.nc','rad_net');
rn_je = nc_varget('je.nc','rad_net');

lwn_jc = nc_varget('jc.nc','lw_net');
lwn_je = nc_varget('je.nc','lw_net');
lwu_jc = nc_varget('jc.nc','lw_up');
lwu_je = nc_varget('je.nc','lw_up');

swn_jc = nc_varget('jc.nc','sw_net');
swn_je = nc_varget('je.nc','sw_net');
swd_jc = nc_varget('jc.nc','sw_down');
swd_je = nc_varget('je.nc','sw_down');
swu_jc = nc_varget('jc.nc','albedo_land');
swu_je = nc_varget('je.nc','albedo_land');


h_jc = nc_varget('jc.nc','ftl_gb');
h_je = nc_varget('je.nc','ftl_gb');
le_jc = nc_varget('jc.nc','latent_heat');
le_je = nc_varget('je.nc','latent_heat');
swc_jc = nc_varget('jc.nc','smcl');
swc_je = nc_varget('je.nc','smcl');
n = 19392;
for i = 1:n
        swc_jc(i,1) = swc_jc(i,1)/100;
        swc_je(i,1) = swc_je(i,1)/100;
        swc_jc(i,2) = swc_jc(i,2)/150;
        swc_je(i,2) = swc_je(i,2)/150;
        swc_jc(i,3) = swc_jc(i,3)/400;
        swc_je(i,3) = swc_je(i,3)/400;
        swc_jc(i,4) = swc_jc(i,4)/1350;
        swc_je(i,4) = swc_je(i,4)/1350;        
end
clear i j n

h_mc = nc_varget('mc.nc','FSH');
h_me = nc_varget('me.nc','FSH');
le_mc = nc_varget('mc.nc','QFX');
le_me = nc_varget('me.nc','QFX');
swc_mc = nc_varget('mc.nc','SMC');
swc_me = nc_varget('me.nc','SMC');

sd_mc = nc_varget('mc.nc','SOLDN');
sd_me = nc_varget('me.nc','SOLDN');
su_mc = nc_varget('mc.nc','FSR');
su_me = nc_varget('me.nc','FSR');
ld_mc = nc_varget('mc.nc','LWDN');
ld_me = nc_varget('me.nc','LWDN');
lu_mc = nc_varget('mc.nc','TRAD');
lu_me = nc_varget('me.nc','TRAD');
n = 38782;
rn_mc = zeros(n,1);
rn_me = zeros(n,1);
for i = 1:n
    rn_mc(i,1) = sd_mc(i,1) + ld_mc(i,1) - su_mc(i,1) - (lu_mc(i,1)^4)*5.67*10^(-8);
    rn_me(i,1) = sd_me(i,1) + ld_me(i,1) - su_me(i,1) - (lu_me(i,1)^4)*5.67*10^(-8);
end
clear i n sd_mc sd_me su_mc su_me ld_mc ld_me lu_mc lu_me



