tic()
%% history
% 2014-01-09 Je-Woo Hong for BVOC  using Guenther et al., 1996


%% information for run
year = 1996; 
    % 1. change the year
    % 2. using ctrl+h 1996 to your year.
    % 3. using ctrl+h control10 to your number.
    % 4. run !

    
%% call sub-material ; for figure and making-grid
lat = nc_varget('control.Gu.nc','latitude');
lon = nc_varget('control.Gu.nc','longitude');
for i = 1:7246
    if lon(1,i)>180
        lon(1,i) = lon(1,i)-360;
    end
end
clear i


%% zenith angle
a = 360*8;
cosz1 = zeros(a,7246);
for i = 1:7246
    lon_z = lon(1,i);
    lambda = lon_z;
    lat_z = lat(1,i);
    fi = lat_z;
    tot_d = a/8;
    temp=0;
    for k = 1:a
        if mod(k,8)==1
            temp = temp + 1;
        end
        tau = 2*pi()*(temp-1)/360;
        delta = 0.006918 -0.399912*cos(tau) +0.070257*sin(tau) -0.006758*cos(2*tau) +0.000907*sin(2*tau);
        Et = 229.18*(0.000075 +0.001868*cos(tau) -0.032077*sin(tau) -0.014615*cos(2*tau) -0.04089*sin(2*tau));
        UT = 90 +(180*(mod(k-1,8)));
        ha = 2*pi*(0.5 +(UT +Et)/1440 +lambda/360);
        cosz1(k,i) = sin(delta)*sin(fi*pi/180)+cos(delta)*cos(fi*pi/180)*cos(ha);
    end
end
clear i k a lon_z lambda lat_z fi tot_d temp tau delta Et UT ha 
cosz = cosz1;
clear cosz1


%% number of days in a month
month_day = [30;30;30;30;30;30;30;30;30;30;30;30];


%% number of 3hr-data in a month ; i.e) SW, T_sfc
month_data3 = zeros(12,1);
for i = 1:12
    month_data3(i,1) = month_day(i,1)*8; % 8 data a day
end
clear i
% at here, took 256 seconds


%% call the LAI and FRAC for this year(1996)
LAI_s = nc_varget('control.Gu.nc','lai'); % monthly-mean lai for each PFT
    % dimension : 12(month), 5(pft), 7246(number of land)
FRAC_s = nc_varget('pft_frac.nc','field1391'); % fraction for eacg tile ; constant
    % dimension : 9(tile), 145, 192
LAI = zeros(12,5,145,192);
for i = 1:12
    for j = 1:240
        for k = 1:7246
            for l = 1:5
                LAI(i,l,add(1,k),add(2,k)) = LAI_s((i-1)*240+j,k);
            end
        end
    end
end
clear i j k l LAI_s
for i = 1:12
    for j = 1:5
        for k = 1:145
            for l = 1:192
                LAI(i,j,k,l) = LAI(i,j,k,l)/240;
            end
        end
    end
end
clear i j k l 

FRAC = zeros(12,5,145,192);
for i = 1:12
    for j = 1:5
        for k = 1:145
            for l = 1:192
                FRAC(i,j,k,l) = FRAC_s(j,k,l);
            end
        end
    end
end
clear i j k l FRAC_s
% only 2.5 seconds for this section


%% main-run
tic()
% for result
isoprene1996 = zeros(13,360,720); % 1-12, and annual
terpene1996 = zeros(13,360,720); % 1-12, and annual

% constant
rgc = 8.314; % molar gas const
alpha = 0.0027;
cli = 1.066;
ct1 = 9.5*10^4; 
ct2 = 2.3*10^5;
ts = 303.0;
tmax = 314.0;
beta = 0.09;
cos_mlsa = 0.5; % Cos(Mean Leaf-Sun Angle); assumed 60 deg.
w2photon = 0.219;
iefac = [24.0; 8.0; 16.0; 16.0; 16.0];
tefac = [0.4; 2.4; 0.8; 0.8; 0.8];
slw = [125.0; 150.0; 125.0; 125.0; 125.0];
red_fac = 0.6875; % Reduce emissions to 550 Tg/yr
% na=6.022e23; % Avogadro no. molecules/mol
% mcarb=0.012;
% icf = na * 1.0e-9 / (5.0 * mcarb * 3600.0);
% tcf = na * 1.0e-9 / (10.0 * mcarb * 3600.0);
% run 

wait = waitbar(0,'...calculation...');
SW_s = nc_varget('control.Gu.nc','sw_down'); % grid downward SW radiation
    % dimension : 3hr, 7246 ; from JULES RUN
T_s = nc_varget('control.Gu.nc','tstar'); % surface temperature for PFT
    % dimension : 3hr + "1", 9, 7246 ; from JULES RUN
    % call SW and surface temperature
SW = zeros(2880, 145, 192);
T = zeros(2880, 5, 145, 192);
for i = 1:2880
    for j = 1:7246
        SW(i,add(1,j),add(2,j)) = SW_s(i,j);
        for k = 1:5
            T(i,k,add(1,j),add(2,j)) = T_s(i,k,j);
        end
    end
end
clear i j k SW_s T_s

for i = 1:12
    
    
    
    temp_data3 = 0;
    for j = 1:7246
        waitbar(((j+(i-1)*7246)/(7246*12)),wait,sprintf('%f',(j+(i-1)*7246)/(7246*12)*100))
        temp =0; % test for no-leaf
        x = add_result(j,1);
        y = add_result(j,2);
        for k = 1:5
            temp = temp + LAI(i,k,x,y);
        end
        clear k
        if temp <=0 % means no-leaf
            % pass !
        else
            if j==1
                temp_data3_2 = temp_data3;
            else
                temp_data3 = temp_data3_2;
            end
            for k = 1:month_data3(i,1)
                temp_data3 = temp_data3 +1;
                TOTPAR = SW(k,x,y)*0.5;
                DIRPAR = TOTPAR*0.7;
                cza = cosz(temp_data3,j);
                pbeam = DIRPAR/w2photon;
                pscat = (TOTPAR - DIRPAR)/w2photon;
                fi = 1;
                ft = 0;
                for l = 1:5 % pft
                    if FRAC(i,l,x,y)>0
                        lw = slw(l,1);
                        ifac = iefac(l,1);
                        tfac = tefac(l,1);
                        tt = T(k,l,x,y);
                        al = LAI(i,l,x,y);
                        if cza>0
                            lai_sun = (1.0 -exp(-0.5*al/cza))*cza/cos_mlsa;
                            lai_shade = al -lai_sun;
                            q1 = 0.07*pbeam*(1.1 -0.1*al)*exp(-cza);
                            qshade = pscat *exp(-0.5*(al^0.7)) +q1;
                            qsun = pbeam*(cos_mlsa/cza) +qshade;
                            v = alpha*qsun;
                            clsun = cli*v/sqrt(1+v^2);
                            v = alpha*qshade;
                            clshade = cli*v/sqrt(1+v^2);
                            denom = rgc*ts*tt;
                            ct = exp((ct1*(tt-ts)/denom))/(1 +exp((ct2*(tt-tmax)/denom)));
                            fi = fi +lw*ct*ifac*FRAC(i,l,x,y)*((lai_sun*clsun) +(lai_shade*clshade));
                            isoprene1996(i,x,y) = isoprene1996(i,x,y)+ fi*red_fac;
                            isoprene1996(13,x,y) = isoprene1996(13,x,y)+ fi*red_fac;
                        end
                        gamma = exp(beta*(tt -ts));
                        ft = ft +al*lw*gamma*tfac*FRAC(i,l,x,y);
                        terpene1996(i,x,y) = terpene1996(i,x,y)+ ft*red_fac;
                        terpene1996(13,x,y) = terpene1996(13,x,y)+ ft*red_fac;
                    end
                end
                clear l
            end
            clear k
        end
    end
    clear j
    
end
close(wait);
clear i
clear na mcarb rgc alpha cli ct1 ct2 ts tmax beta cos_mlsa w2photon red_fac iefac tefac slw icf tcf gamma
clear x y TOTPAR DIRPAR cza pbeam pscat fi ft lw ifac tfac tt al lai_sun lai_shade q1 qshade qsun v clsun clshade denom ct
clear FRAC LAI SW T 
clear temp temp_data3 temp_data3_2 wait
toc()
% for this section ~12800 seconds need !!


