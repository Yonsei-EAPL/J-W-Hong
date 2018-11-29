tic()
%% history
% 2014-01-09 Je-Woo Hong for BVOC  using Guenther et al., 1996

%% information for run
year = 2010; 
    % 1. change the year
    % 2. using ctrl+h 2010 to your year.
    % 3. using ctrl+h control10 to your number.
    % 4. run !

%% call sub-material ; for figure and making-grid
open('submaterial.mat')
m_lat = ans.m_lat; % for figure
m_lon = ans.m_lon; % for figure
add_result = ans.add_result; % for making-grid
area = ans.area; % area for latitude
n_land = ans.n_land; % number of land 
cosz1 = ans.cosz1;
cosz2 = ans.cosz2;
clear ans

%     %% zenith angle
%     a = 365*8;
%     cosz1 = zeros(a,67209);
%     for i = 1:67209
%         lon_z = (add_result(i,2)-1)*0.5 -179.75;
%         lambda = lon_z;
%         lat_z = (add_result(i,1)-1)*0.5 -89.75;
%         fi = lat_z;
%         tot_d = a/8;
%         temp=1;
%         for k = 1:a
%             if mod(k,8)==1
%                 temp = temp + 1;
%             end
%             tau = 2*pi()*(temp-1)/365;
%             delta = 0.006918 -0.399912*cos(tau) +0.070257*sin(tau) -0.006758*cos(2*tau) +0.000907*sin(2*tau);
%             Et = 229.18*(0.000075 +0.001868*cos(tau) -0.032077*sin(tau) -0.014615*cos(2*tau) -0.04089*sin(2*tau));
%             UT = 90 +(180*(mod(k-1,8)));
%             ha = 2*pi*(0.5 +(UT +Et)/1440 +lambda/360);
%             cosz1(k,i) = sin(delta)*sin(fi*pi/180)+cos(delta)*cos(fi*pi/180)*cos(ha);
%         end
%     end
%     clear i k a lon_z lambda lat_z fi tot_d temp tau delta Et UT ha 
%     a = 366*8;
%     cosz2 = zeros(a,67209);
%     for i = 1:67209
%         lon_z = (add_result(i,2)-1)*0.5 -179.75;
%         lambda = lon_z;
%         lat_z = (add_result(i,1)-1)*0.5 -89.75;
%         fi = lat_z;
%         tot_d = a/8;
%         temp=1;
%         for k = 1:a
%             if mod(k,8)==1
%                 temp = temp + 1;
%             end
%             tau = 2*pi()*(temp-1)/366;
%             delta = 0.006918 -0.399912*cos(tau) +0.070257*sin(tau) -0.006758*cos(2*tau) +0.000907*sin(2*tau);
%             Et = 229.18*(0.000075 +0.001868*cos(tau) -0.032077*sin(tau) -0.014615*cos(2*tau) -0.04089*sin(2*tau));
%             UT = 90 +(180*(mod(k-1,8)));
%             ha = 2*pi*(0.5 +(UT +Et)/1440 +lambda/360);
%             cosz2(k,i) = sin(delta)*sin(fi*pi/180)+cos(delta)*cos(fi*pi/180)*cos(ha);
%         end
%     end
%     clear i k a lon_z lambda lat_z fi tot_d temp tau delta Et UT ha 

%% number of days in a month
if mod(year,4)==0
    month_day = [31;29;31;30;31;30;31;31;30;31;30;31];
    cosz = cosz2;
    clear cosz1 cosz2
else
    month_day = [31;28;31;30;31;30;31;31;30;31;30;31];    
    cosz = cosz1;
    clear cosz1 cosz2 
end
clear year 

%% number of 3hr-data in a month ; i.e) SW, T_sfc
month_data3 = zeros(12,1);
for i = 1:12
    month_data3(i,1) = month_day(i,1)*8; % 8 data a day
end
clear i
% at here, took 256 seconds

%% call the LAI and FRAC for this year(2010)
LAI_s = nc_varget('control10.monthly_input.nc','lai'); % monthly-mean lai for each PFT
    % dimension : 12(month), 5(pft), 67209(number of land)
FRAC_s = nc_varget('control10.monthly_input.nc','frac'); % monthly-mean cover fraction for eacg tile
    % dimension : 12, 9(tile), 67209
LAI = zeros(12,5,360,720);
for i = 1:12
    for j = 1:67209
        for k = 1:5
            LAI(i,k,add_result(j,1),add_result(j,2)) = LAI_s(i,k,j);
        end
    end
end
clear i j k LAI_s
FRAC = zeros(12,5,360,720);
for i = 1:12
    for j = 1:67209
        for k = 1:5
            FRAC(i,k,add_result(j,1),add_result(j,2)) = FRAC_s(i,k,j);
        end
    end
end
clear i j k FRAC_s
% only 2.5 seconds for this section

%% main-run
tic()
% for result
isoprene2010 = zeros(13,360,720); % 1-12, and annual
terpene2010 = zeros(13,360,720); % 1-12, and annual

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
for i = 1:12
    % call SW and surface temperature
    if i==1
        SW = nc_varget('SWdown_WFDEI_201001.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_01_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s
    elseif i==2
        SW = nc_varget('SWdown_WFDEI_201002.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_02_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s        
    elseif i==3
        SW = nc_varget('SWdown_WFDEI_201003.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_03_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==4
        SW = nc_varget('SWdown_WFDEI_201004.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_04_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==5
        SW = nc_varget('SWdown_WFDEI_201005.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_05_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==6
        SW = nc_varget('SWdown_WFDEI_201006.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_06_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==7
        SW = nc_varget('SWdown_WFDEI_201007.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_07_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==8
        SW = nc_varget('SWdown_WFDEI_201008.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_08_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==9
        SW = nc_varget('SWdown_WFDEI_201009.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_09_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==10
        SW = nc_varget('SWdown_WFDEI_201010.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_10_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==11
        SW = nc_varget('SWdown_WFDEI_201011.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_11_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    elseif i==12
        SW = nc_varget('SWdown_WFDEI_201012.nc','SWdown'); % grid downward SW radiation
            % dimension : 3hr, 360, 720 ; from WFDEI / JULES community
        T_s = nc_varget('control10.2010_12_3hr.nc','tstar'); % surface temperature for PFT
            % dimension : 3hr + "1", 9, 67209 ; from JULES RUN
        T = zeros(month_data3(i,1), 5, 360, 720);
        for l = 1:month_data3(i,1)
            for j = 1:67209
                for k = 1:5
                    T(l,k,add_result(j,1),add_result(j,2)) = T_s(l,k,j);
                end
            end
        end
        clear l j k T_s                
    end
    temp_data3 = 0;
    for j = 1:67209
        waitbar(((j+(i-1)*67209)/(67209*12)),wait,sprintf('%f',(j+(i-1)*67209)/(67209*12)*100))
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
                            isoprene2010(i,x,y) = isoprene2010(i,x,y)+ fi*red_fac;
                            isoprene2010(13,x,y) = isoprene2010(13,x,y)+ fi*red_fac;
                        end
                        gamma = exp(beta*(tt -ts));
                        ft = ft +al*lw*gamma*tfac*FRAC(i,l,x,y);
                        terpene2010(i,x,y) = terpene2010(i,x,y)+ ft*red_fac;
                        terpene2010(13,x,y) = terpene2010(13,x,y)+ ft*red_fac;
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


