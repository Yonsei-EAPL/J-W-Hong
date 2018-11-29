tic()
%% history
% 2014-01-09 Je-Woo Hong for BVOC  using Guenther et al., 1996
% 2014-08  Woomi Jung : modify for HadGEM data (grid and number of date)

%% information for run
year = 1997; 
    % 1. change the year and filename
    % 2. using ctrl+h 2002 to your year.
    % 3. using ctrl+h control2 to your number.
    % 4. run !
length = 7246; % jules output(using AMIP data) demension length
%% call sub-material ; for figure and making-grid
open('BVOC_info.mat');
add = ans.add; % for making-grid
area = ans.area; % area for latitude
land = ans.land; % number of land 
% cosz1 = ans.cosz1; % 3hr zenith angle (8*360 = 2880)
% cosz3 = ans.cosz3; % 1hr zenith angle (24*360 = 8640)
clear ans

%%

filename = 'F:\BVOC\1997\control.Gu.nc';

lai = zeros(145,192);
temp = zeros(145,192);
for i = 1:145
    for j = 1:192
        lai(i,j) = lai_t(100,2,i,j);
        temp(i,j) = t_t(100,2,i,j);
    end
end
clear i j

%% reference lat/lon calculation for lat/lon degree -> grid point 

lat = nc_varget(filename,'latitude');
lon = nc_varget(filename,'longitude');

% lat_ref = zeros(145,1);
% lat_ref(1,1) = -90;
% for i=1:144
%     lat_ref(i+1,1) = lat_ref(i,1) + 1.25;
% end
% clear i
% 
% lon_ref = zeros(192,1);
% lon_ref(1,1) = 0;
% for i=1:191
%     lon_ref(i+1,1) = lon_ref(i,1) + 1.875;
% end
% clear i

% add = zeros(2,length);
% for i = 1:length
%     add(1,i) = find(lat(i,1)==lat_ref(:,1));
%     add(2,i) = find(lon(i,1)==lon_ref(:,1));
% end
% clear i
 
% area = zeros(145,1);
% for i = 1:145
%     area(i,1) = 138000*6370000*cos(pi()/180*lat_ref(i,1))*2*pi()/360*1.875;
% end
% clear i 

%% zenith angle

% 360 days(HadGEM) hourly zenith angle

     a = 360*24; % 360 day * 24hr
     cosz3 = zeros(a,length);
     for i = 1:length
         lon_z = lon(i,1) -179.0625;
         lambda = lon_z;
         lat_z = lat(i,1);
         fi = lat_z;
         tot_d = a/24;
         temp=0;
         for k = 1:a
             if mod(k,24)==1
                 temp = temp + 1;
             end
             tau = 2*pi()*(temp-1)/360;
             delta = 0.006918 -0.399912*cos(tau) +0.070257*sin(tau) -0.006758*cos(2*tau) +0.000907*sin(2*tau);
             Et = 229.18*(0.000075 +0.001868*cos(tau) -0.032077*sin(tau) -0.014615*cos(2*tau) -0.04089*sin(2*tau));
             UT = 30 +(60*(mod(k-1,24)));
             ha = 2*pi*(0.5 +(UT +Et)/1440 +lambda/360);
             cosz3(k,i) = sin(delta)*sin(fi*pi/180)+cos(delta)*cos(fi*pi/180)*cos(ha);
         end
     end
     clear i k a lon_z lambda lat_z fi tot_d temp tau delta Et UT ha lon lat

% 366 days, 3hr zenith angle

%      a = 366*8;
%      cosz2 = zeros(a,67209);
%      for i = 1:67209
%          lon_z = (add_result(i,2)-1)*0.5 -179.75;
%          lambda = lon_z;
%          lat_z = (add_result(i,1)-1)*0.5 -89.75;
%          fi = lat_z;
%          tot_d = a/8;
%          temp=1;
%          for k = 1:a
%              if mod(k,8)==1
%                  temp = temp + 1;
%              end
%              tau = 2*pi()*(temp-1)/366;
%              delta = 0.006918 -0.399912*cos(tau) +0.070257*sin(tau) -0.006758*cos(2*tau) +0.000907*sin(2*tau);
%              Et = 229.18*(0.000075 +0.001868*cos(tau) -0.032077*sin(tau) -0.014615*cos(2*tau) -0.04089*sin(2*tau));
%              UT = 90 +(180*(mod(k-1,8)));
%              ha = 2*pi*(0.5 +(UT +Et)/1440 +lambda/360);
%              cosz2(k,i) = sin(delta)*sin(fi*pi/180)+cos(delta)*cos(fi*pi/180)*cos(ha);
%          end
%      end
%      clear i k a lon_z lambda lat_z fi tot_d temp tau delta Et UT ha 

% hourly zenith angle -> 3hr mean
temp = 0;
cosz1 = zeros(2880,length);
for i = 1:length
    t = 0;
    for j = 1:8640
        temp = temp+cosz3(j,i);
        if mod(j,3)==0
            t = t+1;
            cosz1(t,i) = temp/3;
            temp = 0;
        end
    end
end
clear i j t temp

% test

cosz4 = zeros(145,192);
for i = 1:7246
    cosz4(add(1,i),add(2,i)) = cosz(3,i);
end
clear i


%% number of days in a month
% if mod(year,4)==0
%     month_day = [31;29;31;30;31;30;31;31;30;31;30;31];
%     cosz = cosz2;
%     clear cosz1 cosz2
% else
%     month_day = [31;28;31;30;31;30;31;31;30;31;30;31];    
%     cosz = cosz1;
%     clear cosz1 cosz2 
% end
% clear year 

gm_1 = zeros(145,192);
for i = 1:7246
    gm_1(add(1,i),add(2,i)) = gm(6,180,i);
end
clear i

% AMIP data : 360 day (30 day * 12 month)

month_day = [30;30;30;30;30;30;30;30;30;30;30;30];
cosz = cosz3;
clear cosz3 year

%% number of 3hr-data in a month ; i.e) SW, T_sfc
month_data3 = zeros(12,1);
for i = 1:12
    month_data3(i,1) = month_day(i,1)*8; % 8 data a day
end
clear i

% yrs = zeros(10,1); % total period(year) of data set - historical : 1996-2005 (10yrs) / RCP : 2090-2099 (10yrs)
% yrs(1,1) = year; % start year -1
% for i=2:10
%     yrs(i,1) = yrs(i-1,1) + 1;
% end
% clear i

%% call the LAI and FRAC for this year(2002)
LAI_s = nc_varget(filename,'lai'); % monthly-mean lai for each PFT
    % dimension : 12(month)*10yr, 5(pft), 1, 7246(number of land)
FRAC_s = nc_varget('pft_frac.nc','field1391'); % monthly-mean cover fraction for eacg tile (constant)
    % dimension : 9(tile), 145*192

% read variables as 2-D [145 192]
% lai from 1996-2005 monthly data set -> nt = 120 (1996:1~12, 1997: 13~24)

% LAI = zeros(12,5,145,192);
% for i = 1:12
%     for j = 1:length
%         for k = 1:5
%             LAI(i,k,find(lat(j,1)==lat_ref(:,1)),find(lon(j,1)==lon_ref(:,1))) = LAI_s(i+12*find(year==yrs(:))-1,k,1,j);
%         end
%     end
% end
% clear i j k LAI_s

% % LAI 3hr data -> monthly data (monthly mean) : 3hr-data 8*30(day) = 240
LAI = zeros(12,5,145,192);
for i = 1:12
    for l = 1:5
        for k = 1:length
            for j = 1:240
                LAI(i,l,add(1,k),add(2,k)) = LAI(i,l,add(1,k),add(2,k))+ LAI_s((i-1)*240+j,l,k);
            end
        end
    end
end
clear i j k l
for i = 1:12
    for j = 1:145
        for k = 1:192
            for l = 1:5
                LAI(i,l,j,k) = LAI(i,l,j,k)/240;
            end
        end
    end
end
clear i j k l LAI_s

% lai test

frac_t = zeros(145,192);
for i = 1:145
    for j = 1:192
        frac_t(i,j) = FRAC(1,5,i,j);
    end
end
clear i j

% read fraction data

FRAC = zeros(12,5,145,192);
for i = 1:12
    for j = 1:145
        for t = 1:192
            for k = 1:5
                FRAC(i,k,j,t) = FRAC_s(k,j,t);
            end
        end
    end
end
clear i j k t FRAC_s

toc()
%% main-run
tic()
% for result : change the year and shift+enter
isoprene1997 = zeros(13,145,192); % 1-12, and annual
terpene1997 = zeros(13,145,192);
gm = zeros(12,240,length);
% 1-12, and annual

% constant
rgc = 8.314; % molar gas const [J mol-1 K-1]
alpha = 0.0027;
cli = 1.066;
ct1 = 9.5*10^4; % [J mol-1]
ct2 = 2.3*10^5; % [J mol-1]
ts = 303.0; % [K]
tmax = 314.0; % [K]
beta = 0.09; % 0.057~0.144 -> 0.09 +- 0.015[K-1]
cos_mlsa = 0.5; % Cos(Mean Leaf-Sun Angle); assumed 60 deg.
w2photon = 0.219; % = na * h(Plank's const.) * c / lambda [J umol-1]
iefac = [24.0; 8.0; 16.0; 16.0; 16.0];
tefac = [0.4; 2.4; 0.8; 0.8; 0.8];
slw = [125.0; 150.0; 125.0; 125.0; 125.0];
red_fac = 0.6875; % Reduce emissions to 550 Tg/yr
% na=6.022e23; % Avogadro no. molecules/mol
% mcarb=0.012;
% icf = na * 1.0e-9 / (5.0 * mcarb * 3600.0);
% tcf = na * 1.0e-9 / (10.0 * mcarb * 3600.0);

SW_s = nc_varget(filename,'sw_down'); % SW radiation
% dimension : 3hr*2880 + 1, 7246 ; from JULES RUN
T_s = nc_varget(filename,'tstar'); % surface temperature for PFT
% dimension : 3hr*2880 + 1, 9, 7246 ; from JULES RUN


% SW = zeros(2880,145,192);
% for i=1:12
%     lst=(i-1)*month_data3(i,1);
%     
% for l = 1:2880
%     for j = 1:length
%         SW(l,add(1,j),add(2,j)) = SW_s(l,1,j);
%     end
% end
% clear l j
% 
% tic()
% run 
for i = 1:12
    lst = (i-1)*month_data3(i,1);
    % call SW and surface temperature        
    T = zeros(month_data3(i,1), 5, 145, 192);
    n = 0;
    for l = lst+1:lst+month_data3(i,1)
        n = n+1;
        for j = 1:length
            for k = 1:5
                T(n,k,add(1,j),add(2,j)) = T_s(l,k,1,j);
            end
        end
    end
    clear l j k n
    SW = zeros(month_data3(i,1),145,192);
    n = 0;
    for l = lst+1:lst+month_data3(i,1)
        n = n+1;
        for j = 1:length
            SW(n,add(1,j),add(2,j)) = SW_s(l,1,j);
        end
    end
    clear l j n
  
%    temp_data3 = 0;
    for j = 1:length
        temp =0; % test for no-leaf
        x = add(1,j);
        y = add(2,j);
        for k = 1:5
            temp = temp + LAI(i,k,x,y);
        end
        clear k
        if temp <=0 % means no-leaf
            % pass !
        else
            for k = 1:month_data3(i,1)
                TOTPAR = SW(k,x,y)*0.5;
                DIRPAR = TOTPAR*0.7;                
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
                        for t = 1:3
                            cza = cosz((k-1)*3+t,j);
                            if cza>0
                                lai_sun = (1.0 - exp((-0.5)*al/cza))*cza/cos_mlsa;
                                lai_shade = al - lai_sun;
                                q1 = 0.07*pbeam*(1.1 - 0.1*al)*exp((-1.0)*cza);
                                qshade = pscat*exp((-0.5)*(al^(0.7))) + q1;
                                qsun = pbeam*(cos_mlsa/cza) + qshade;
                                v = alpha*qsun;
                                clsun = cli*v/sqrt(1.0+v^2);
                                v = alpha*qshade;
                                clshade = cli*v/sqrt(1.0+v^2);
                                denom = rgc*ts*tt;
                                ct = (exp(ct1*(tt-ts)/denom))/(1.0+exp(ct2*(tt-tmax)/denom));                                
                                fi = fi + FRAC(i,l,x,y)*lw*ifac*ct*((lai_sun*clsun)+(lai_shade*clshade));
                                if l==1
                                    % isoprene gamma test
                                    gm(i,k,j) = ct*((lai_sun*clsun)+(lai_shade*clshade));
                                end
                            end
                        end
                        clear t
                        gamma = exp(beta*(tt -ts));
                        ft = ft +al*lw*gamma*tfac*FRAC(i,l,x,y);
                    end                
                end
                isoprene1997(i,x,y) = isoprene1997(i,x,y)+ fi*red_fac;
                isoprene1997(13,x,y) = isoprene1997(13,x,y)+ fi*red_fac;
                terpene1997(i,x,y) = terpene1997(i,x,y)+ ft*red_fac;
                terpene1997(13,x,y) = terpene1997(13,x,y)+ ft*red_fac;
                clear l
            end
            clear k n
        end
    end
    clear j
end
clear i 
clear na mcarb rgc alpha cli ct1 ct2 ts tmax beta cos_mlsa w2photon red_fac iefac tefac slw icf tcf 
clear x y TOTPAR DIRPAR cza pbeam pscat fi ft lw ifac tfac tt al lai_sun lai_shade q1 qshade qsun v clsun clshade denom ct
clear LAI SW T SW_s T_s
toc()
% for this section ~2300 seconds (~40 minutes) need !!

% open('결과정리.mat');
% area_jw = ans.area;
% clear ans

% annual data read & unit convert(ug/m2/hr -> *10^(-6)(g/ug)*area(m2)*3(hr))

iso = zeros(145,192);
ter = zeros(145,192);
for i = 1:145
    for j = 1:192
        iso(i,j) = isoprene1997(13,i,j)*area(i)/1000000;
        ter(i,j) = terpene1997(13,i,j)*3*area(i)/1000000;
        if ter(i,j)>10^11 % 
            ter(i,j) = 0;
        end
    end
end
clear i j 

% T test
ts = zeros(145,192);
for i = 1:145
    for j = 1:192
        ts(i,j) = T(5,4,i,j);
    end
end
clear i j

% the amount of annual (total) emission
iso(~isfinite(iso)) = 0;
iso1996 = 0;
for i = 1:145
    for j = 1:192
        iso1996 = iso1996 + iso(i,j);
    end
end
clear i j

% calculation of hist/rcp mean

iso_hist = zeros(10,145,192);
iso_rcp = zeros(10,145,192);
ter_hist = zeros(10,145,192);
ter_rcp = zeros(10,145,192);

for i = 1:145
    for j = 1:192
        iso_hist(1,i,j) = isoprene1997(13,i,j);
        iso_hist(2,i,j) = isoprene1997(13,i,j);
        iso_hist(3,i,j) = isoprene1998(13,i,j);
        iso_hist(4,i,j) = isoprene1999(13,i,j);
        iso_hist(5,i,j) = isoprene2000(13,i,j);
        iso_hist(6,i,j) = isoprene2001(13,i,j);
        iso_hist(7,i,j) = isoprene2002(13,i,j);
        iso_hist(8,i,j) = isoprene2003(13,i,j);
        iso_hist(9,i,j) = isoprene2004(13,i,j);
        iso_hist(10,i,j) = isoprene2005(13,i,j);
        
        iso_rcp(1,i,j) = isoprene2090(13,i,j);
        iso_rcp(2,i,j) = isoprene2091(13,i,j);
        iso_rcp(3,i,j) = isoprene2092(13,i,j);
        iso_rcp(4,i,j) = isoprene2093(13,i,j);
        iso_rcp(5,i,j) = isoprene2094(13,i,j);
        iso_rcp(6,i,j) = isoprene2095(13,i,j);
        iso_rcp(7,i,j) = isoprene2096(13,i,j);
        iso_rcp(8,i,j) = isoprene2097(13,i,j);
        iso_rcp(9,i,j) = isoprene2098(13,i,j);
        iso_rcp(10,i,j) = isoprene1997(13,i,j);
        
        ter_hist(1,i,j) = terpene1997(13,i,j);
        ter_hist(2,i,j) = terpene1997(13,i,j);
        ter_hist(3,i,j) = terpene1998(13,i,j);
        ter_hist(4,i,j) = terpene1999(13,i,j);
        ter_hist(5,i,j) = terpene2000(13,i,j);
        ter_hist(6,i,j) = terpene2001(13,i,j);
        ter_hist(7,i,j) = terpene2002(13,i,j);
        ter_hist(8,i,j) = terpene2003(13,i,j);
        ter_hist(9,i,j) = terpene2004(13,i,j);
        ter_hist(10,i,j) = terpene2005(13,i,j);
        
        ter_rcp(1,i,j) = terpene2090(13,i,j);
        ter_rcp(2,i,j) = terpene2091(13,i,j);
        ter_rcp(3,i,j) = terpene2092(13,i,j);
        ter_rcp(4,i,j) = terpene2093(13,i,j);
        ter_rcp(5,i,j) = terpene2094(13,i,j);
        ter_rcp(6,i,j) = terpene2095(13,i,j);
        ter_rcp(7,i,j) = terpene2096(13,i,j);
        ter_rcp(8,i,j) = terpene2097(13,i,j);
        ter_rcp(9,i,j) = terpene2098(13,i,j);
        ter_rcp(10,i,j) = terpene1997(13,i,j);
    end
end

clear i j

% calculation of hist/rcp (10yrs) mean
iso_hist_m = zeros(145,192);
iso_rcp_m = zeros(145,192);
ter_hist_m = zeros(145,192);
ter_rcp_m = zeros(145,192);

for i = 1:145
    for j = 1:192
        for t = 1:10
            iso_hist_m(i,j) = iso_hist_m(i,j) + iso_hist(t,i,j);
            iso_rcp_m(i,j) = iso_rcp_m(i,j) + iso_rcp(t,i,j);

            temp1 = ter_hist(t,i,j)*3*area_jw(i)/1000000;
            temp2 = ter_rcp(t,i,j)*3*area_jw(i)/1000000;
            if temp1 > 10^11
                temp1 = 0;
            end
            if temp2 > 10^11
                temp2 = 0;
            end
            ter_hist_m(i,j) = ter_hist_m(i,j) + temp1;
            ter_rcp_m(i,j) = ter_rcp_m(i,j) + temp2;
            
            if t == 10
                iso_hist_m(i,j) = iso_hist_m(i,j)*3*area_jw(i)/10000000;
                iso_rcp_m(i,j) = iso_rcp_m(i,j)*3*area_jw(i)/10000000;
                ter_hist_m(i,j) = ter_hist_m(i,j)/10;
                ter_rcp_m(i,j) = ter_rcp_m(i,j)/10;
            end            
        end
    end
end

clear i j t temp1 temp2

% calculation of difference between hist and rcp

iso_diff = zeros(145,192);
ter_diff = zeros(145,192);

for i = 1:145
    for j = 1:192
        iso_diff(i,j) = iso_rcp_m(i,j) - iso_hist_m(i,j);
        ter_diff(i,j) = ter_rcp_m(i,j) - ter_hist_m(i,j);
    end
end
clear i j