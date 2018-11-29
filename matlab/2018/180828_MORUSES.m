%% input
wrr = (0.61)*0.9; % EP ctl
hwr = 0.78; % EP ctl
hgt = 21.4; % EP ctl

%% constant
z1_uv = 10;
z1_tq = 10;
a = 4.43;
cdz =1.2;
kappa2 = 0.16;
z0m_mat = 0.05;
alpha = 0.15; % get_us
u1 = 1; % get_us
v1 = 0; % get_us
vkman = 0.4; %von Karman constant

%% init_urban
sc_hwr = 0.5*(hwr/ (2 * atan(1)));
d_h = 1.0 - wrr*(a^(wrr-1));
disp = d_h *hgt;
ztm = (1-d_h)*exp(-(cdz*(1-d_h)*sc_hwr*wrr/kappa2)^(-0.5))*hgt;
ztm = max(ztm, z0m_mat);

%% urbanz0
z0m = z0m_mat;

if (z0m-10^(-6) <= ztm)
    zref = 0.1*hgt;
    u_us = log(z1_uv/ztm+1)/sqrt(kappa2);
    
    % for urban-roof
    roof_log = max(1.1*hgt- disp, ztm);
    fz1_roof = log(roof_log/ztm)/ log(z1_uv/ztm +1);
    ures1_roof = log(zref/z0m +1)*log((zref+z0m)/(0.1*z0m));
    ures1_roof = ures1_roof/(kappa2 * fz1_roof);
    ures3_roof = (1- fz1_roof)*(u_us^2);
    urest_roof = ures1_roof + ures3_roof;
    zth_h = kappa2*urest_roof/log(z1_uv/ztm +1);
    zth_h = (z1_tq + ztm)/(exp(zth_h) * hgt);
    zth_h_roof = max(zth_h * hgt, 1e-30);
	% for urban canyon
    % CALL get_us for urest_can
        bht = hgt;
        z0t = ztm;
        z1_u = z1_uv;
        lrw  = 3*(2*hwr/pi());
        ltw  = 1.5*(2*hwr/pi());
        hrh  = (lrw-1)/(lrw-ltw);
        lrw  = min(lrw,1);
        ltw  = min(ltw,1);
        hrh  = max(hrh,0);
        hrh  = max(hrh,1);
        lttw = lrw + (1- hrh) *sqrt((lrw- ltw)^2 + hwr^2);
        ut  = sqrt(u1^2 + v1^2) * 2/pi();
        vt  = sqrt(u1^2 + v1^2) * 2/pi();
        uu1 = sqrt(u1^2 + v1^2);
        log_inf = max(bht-disp, z0t);
        uct = ut * log(log_inf/z0t) / log(z1_u/z0t + 1);
        vct = vt * log(log_inf/z0t) / log(z1_u/z0t + 1);
        alocan1 = vct * log(zref/z0m)/log(bht/z0m);  %street
        alocan2 = vct * (1/(1-(z0m/bht)) - 1/log(bht/z0m)); %wall
        wd = bht/hwr;
        lr = 3*bht;
%         if hwr>2/3
        lt = 0;
        hs = 0;
        
        nu = lr/(2*wd);
        udw = uct*(1-exp(-alpha))/alpha;
        ud = udw;
        uu = uct*bht*exp(-nu*alpha)*(1-exp(-nu*alpha*wd/bht));
        uu = uu/(nu*alpha*wd);
        uuw = uct*exp(-nu*alpha*(bht+wd)/bht)*(1-exp(-nu*alpha));
        uuw = uuw/(nu*alpha);
        uu = sqrt(uu^2 + alocan1^2);	%street
        ud = sqrt(ud^2 + alocan1^2);
        uuw = sqrt(uuw^2 + alocan2^2);	%wall
        udw = sqrt(udw^2 + alocan2^2);
        r5 = (1-uu)*u_us^2/(1+lttw-ltw) /(min(lr/2,wd)/wd);
        z0h = 0.1*z0m;
        r8 = log(zref/z0m)*log(zref/z0h);
        r3 = r8/(uuw* kappa2)/hwr;
        r4 = r8/(uu*kappa2) /(min(lr, wd)/wd);
        r6 = r8/(udw*kappa2);
        r8 = r8/(ud*kappa2);
        r7 = (1-udw)*u_us^2;

        scaler6hat = max(1e-3, (bht-hs)/wd);
        scaler6 = max(1e-3, hs/wd);
        scaler7 = max(1e-3, (wd-min(lr/2, wd))/wd);
        scaler8 = max(1e-3, (wd-min(lr,wd))/wd);

        r6hat = r6/scaler6hat;
        r6 = r6/scaler6;
        r7 = r7/scaler7;
        r8 = r8/scaler8;

        rbulk1 = (r3*r4*r6hat)/((r3*r4)+(r3*r6hat)+(r4*r6hat)) + r5;
        rbulk2 = (r6*r8)/(r6+r8) + r7;

        rbulk = 1/ ((1/rbulk1)+(1/rbulk2));
        urest_can = rbulk;
%         end
    % return to urbanz0
    zth_h = kappa2*urest_can/log(z1_uv/ztm +1);
    zth_h = (z1_tq + ztm)/(exp(zth_h) * hgt);
    zth_h_can = max(zth_h * hgt, 1e-30);
end

z0m_tile = ztm;
z0h_tile_r = zth_h_roof;
z0h_tile_c = zth_h_can;

clear alocan1 alocan2 bht fz1_roof hrh hs log_inf lr lrw lt lttw ltw nu r3 r4 r5 r6 r6hat r7 r8
clear rbulk rbulk1 rbulk2 roof_log scaler6 scaler6hat scaler7 scaler8 u1 u_us
clear uct ud udw ures1_roof ures3_roof ut uu uu1 uuw v1 vct vt wd zref zth_h zth_h_can zth_h_roof


%% sf_exch
wind_profile_factor = 1;

%line 1111
zetam_r = log((z1_uv + z0m_tile)/z0m_tile);
zetam_c = log((z1_uv + z0m_tile)/z0m_tile);
zetah_r = log((z1_tq + z0m_tile)/z0h_tile_r);
zetah_c = log((z1_tq + z0m_tile)/z0h_tile_c);
chn_r = (vkman/zetah_r)*(vkman/zetam_r)*wind_profile_factor;
chn_c = (vkman/zetah_c)*(vkman/zetam_c)*wind_profile_factor;





