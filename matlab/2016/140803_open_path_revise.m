tic
%% History
% 130714 Mr.Je-Woo Hong & Prof.Jinkyu Hong, discussed with Mr.Keunmin Lee 
% 131001 Je-Woo Hong ; for molar volume, using observed pressure/
%                      for open-path data process (using density)
% 140713 Mr.Keunmin Lee; modification to load data continuously

%% load input data
dataDir = 'E:\EAPL\KBSI\es'; % folder name
dataName = dir(fullfile(dataDir, 'CSV_5992.ts_data_*'))'; %for ochang
% dataName = dir(fullfile(dataDir, 'CSV_7679.ts_data_*'))'; %for hongcheon
%dataName = dir(fullfile(dataDir, 'CSV_7681.ts_data_*'))'; %for samcheok
%dataName = dir(fullfile(dataDir, 'CSV_7682.ts_data_*'))'; %for pyeongchang
%dataName = dir(fullfile(dataDir, 'CSV_6330.ts_data_*'))'; %for jeju

total_result=cell(length(dataName),2);

for n=1:length(dataName)
    data=importdata(fullfile(dataDir, dataName(n).name));    
    
%% input(*.csv) information
sonic_ang = 0; %for ochang
% sonic_ang = 90+8.15; %for hongcheon
%sonic_ang = 20+8.13; %for samcheok
%sonic_ang = 0+8.16; %for pyeongchang
%sonic_ang = 0+6.7; %for jeju

[size_n size_var] = size(data);
num_30min = 0;
num_30min_n = zeros(1,3);

for i = 1:size_n
    if i==1
        num_30min = 1;
        num_30min_n(num_30min,1) = num_30min;
        num_30min_n(num_30min,2) = 1;
        num_30min_n(num_30min,3) = 1;
    else
%        if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==0.05)          
        if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==0.1)
%          if ((mod(data(i,3),100)==0)||(mod(data(i,3),10)==0))&&(mod(data(i,4),1)==0.1)
            num_30min = num_30min+1;
            num_30min_n(num_30min,1) = num_30min;
            num_30min_n(num_30min,2) = 1;
            num_30min_n(num_30min,3) = num_30min_n(num_30min-1,3)+1;
        else
            num_30min_n(num_30min,2) = num_30min_n(num_30min,2)+1;
            num_30min_n(num_30min,3) = num_30min_n(num_30min,3)+1;
        end
    end
end
clear i num_30min


%% result

result = zeros(max(num_30min_n(:,1)),54);
% 1; n_data (unitless)
% 2; mean wind speed, U (m/s)
% 3; u_bar (m/s)
% 4; v_bar (m/s)
% 5; w_bar (m/s)
% 6; wind direction (including sonic_angle) (degree)
% 7; std_wd (degree)
% 8; u* (m/s)
% 9; cov_uw (m2/s2)
% 10; alpha (1st, rotation)
% 11; beta (2nd, rotation)
% 12; mean Ts (degree C)
% 13; mean Tk (degree C)
% 14; mean ambient temperature (degree C)
% 15; mean ambient pressure (kPa)
% 16; mean CO2 (umol/mol)
% 17; min CO2 (umol/mol)
% 18; max CO2 (umol/mol)
% 19; mean H2O (mmol/mol)
% 20; min H2O (mmol/mol)
% 21; max H2O (mmol/mol)
% 22; e (kPa)
% 23; es (kPa)
% 24; RH (%)
% 25; rho_d (g/m3)
% 26; rho_v (gH2O/m3)
% 27; rho_co2 (mgCO2/m3)
% 28; rhoCp (J/K/m3)
% 29; cov_wTs (K*m/s)
% 30; rhoCp*cov_wTs (W/m2)
% 31; Qh ; AmeriFlux (W/m2) correct the effect of h2o on Ts
% 32; lambda ; latent heat of vaporization of water (J/g)
% 33; Qe (W/m2)
% 34; Fc (umol/m2/s)
% 35; std_u (m/s)
% 36; std_v (m/s)
% 37; std_w (m/s)
% 38; std_u/u* (unitless)
% 39; std_v/u* (unitless)
% 40; std_w/u* (unitless)
% 41; std_u/U (unitless)
% 42; std_v/U (unitless)
% 43; std_w/U (unitless)
% 44; r_uw (unitless)
% 45; sk_u (unitless)
% 46; sk_v (unitless)
% 47; sk_w (unitless)
% 48; sk_Ts (unitless)
% 49; ku_u (unitless)
% 50; ku_v (unitless)
% 51; ku_w (unitless)
% 52; ku_Ts (unitless)
% 53; Qe_WPL (W/m2)
% 54; Fc_WPL (umol/m2/s)


%% position

po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8;
po_CO2 = 10;
po_H2O = 11;
po_amb_tmp = 13;
po_amb_prs = 14;


%% constant
R = 8.3143*10^(-6); % kPa*m3/K/mmol
%Mc = 12; % mgC/mmol
Mco2 = 0.04401; % mg/umol
Md = 0.029; % g/mmol
Mv = 0.018015; % g/mmol
    

%% main process
for i = 1:max(num_30min_n(:,1))
    result(i,1) = num_30min_n(i,2);    
    
    % extract
    temp = zeros(num_30min_n(i,2),8);
    for j = 1:num_30min_n(i,2)
        if i == 1
            temp(j,1) = data(j,po_u);
            temp(j,2) = data(j,po_v);
            temp(j,3) = data(j,po_w);
            temp(j,4) = data(j,po_Ts);
            temp(j,5) = data(j,po_CO2);
            temp(j,6) = data(j,po_H2O);
            temp(j,7) = data(j,po_amb_prs);
            temp(j,8) = data(j,po_amb_tmp);
        else
            temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
            temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
            temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
            temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
            temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
            temp(j,6) = data(num_30min_n(i-1,3)+j,po_H2O);
            temp(j,7) = data(num_30min_n(i-1,3)+j,po_amb_prs);
            temp(j,8) = data(num_30min_n(i-1,3)+j,po_amb_tmp);
        end
    end
    clear j
    
    % mean wind-speed, U
    temp_ws = 0;
    for j = 1:num_30min_n(i,2)
        temp_ws = temp_ws + (temp(j,1)^2 + temp(j,2)^2 + temp(j,3)^2)^(0.5);
    end
    temp_ws = temp_ws/num_30min_n(i,2);
    result(i,2) = temp_ws; % result : mean wind-speed, U
    clear j temp_ws
    
    % mean wind-direction
    temp_std_wd = 0;
    u_bar = mean(temp(:,1));
    v_bar = mean(temp(:,2));
    w_bar = mean(temp(:,3));
    if u_bar>0
        if v_bar>0
            temp_wd = 360 - atan(v_bar/u_bar)/pi()*180;
        else
            temp_wd = atan((-1*v_bar)/u_bar)/pi()*180;
        end
    else
        if v_bar>0
            temp_wd = 180 + atan(v_bar/(-1*u_bar))/pi()*180;
        else
            temp_wd = 180 - atan(v_bar/u_bar)/pi()*180;
        end
    end
    mean_wd = temp_wd;
    for j = 1:num_30min_n(i,2)
        if temp(j,1)>0
            if temp(j,2)>0
                temp_wd = 360 - atan(temp(j,2)/temp(j,1))/pi()*180;
            else
                temp_wd = atan((-1*temp(j,2))/temp(j,1))/pi()*180;
            end
        else
            if temp(j,2)>0
                temp_wd = 180 + atan(temp(j,2)/(-1*temp(j,1)))/pi()*180;
            else
                temp_wd = 180 - atan(temp(j,2)/temp(j,1))/pi()*180;
            end
        end
        if abs(mean_wd-temp_wd)>180
            if (mean_wd-temp_wd)<0
                temp_std_wd = temp_std_wd + (mean_wd-temp_wd+360)^2;
            else
                temp_std_wd = temp_std_wd + (360 - mean_wd-temp_wd)^2;
            end
        else
            temp_std_wd = temp_std_wd + (mean_wd-temp_wd)^2;
        end
    end
    temp_std_wd = (temp_std_wd/num_30min_n(i,2))^(0.5);
    result(i,7) = temp_std_wd; % result : standard deviation of wind-direction
    if mean_wd+sonic_ang>360
        mean_wd = mean_wd + sonic_ang-360;
    else
        mean_wd = mean_wd + sonic_ang;
    end
    result(i,6) = mean_wd; % result : wind direction
    clear j temp_wd temp_std_wd mean_wd 

    % double rotation
    temp_wind = zeros(num_30min_n(i,2),3);
    alpha = atan2(v_bar,u_bar);
    beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.5)));
    for j = 1:num_30min_n(i,2)
        temp_wind(j,1) = cos(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+sin(beta)*temp(j,3);
        temp_wind(j,2) = -sin(alpha)*temp(j,1)+cos(alpha)*temp(j,2);
        temp_wind(j,3) = -sin(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+cos(beta)*temp(j,3);
        temp(j,1) = temp_wind(j,1);
        temp(j,2) = temp_wind(j,2);
        temp(j,3) = temp_wind(j,3);
    end
    result(i,10)= alpha; % result : 1st rotation angle
    result(i,11)= beta; % result : 2nd rotation angle
    result(i,35)=std(temp(:,1)); % result : standard deviation of u
    result(i,36)=std(temp(:,2)); % result : standard deviation of v
    result(i,37)=std(temp(:,3)); % result : standard deviation of w
   result(i,45)=skewness(temp(:,1)); % result : skewness of u
   result(i,49)=kurtosis(temp(:,1)); % result : kurtosis of u
   result(i,46)=skewness(temp(:,2)); % result : skewness of v
   result(i,50)=kurtosis(temp(:,2)); % result : kurtosis of v
   result(i,47)=skewness(temp(:,3)); % result : skewness of w
   result(i,51)=kurtosis(temp(:,3)); % result : kurtosis of w    
    clear u_bar v_bar w_bar alpha beta j temp_wind

    % friction velocity
        % u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.25) 
        % Stull, p.67
    uw = 0; % for mean(u'w')
    vw = 0; % for mean(v'w')
    u_bar = mean(temp(:,1));
    result(i,3) = u_bar; % result : mean u
    v_bar = mean(temp(:,2));
    result(i,4) = v_bar; % result : mean v
    w_bar = mean(temp(:,3));
    result(i,5) = w_bar; % result : mean w
    for j = 1:num_30min_n(i,2)
        uw = uw + (temp(j,1)-u_bar)*(temp(j,3)-w_bar);
        vw = vw + (temp(j,2)-v_bar)*(temp(j,3)-w_bar);
    end
    uw = uw/num_30min_n(i,2);
    vw = vw/num_30min_n(i,2);
    result(i,8) = (uw^2 + vw^2)^(0.25); % result : u_star
    result(i,9) = uw; % result : cov_uw
    result(i,44) = uw/result(i,35)/result(i,37);
    clear uw vw u_bar v_bar w_bar j
    
    % std_u,v,w/u*,U
    result(i,38) = result(i,35)/result(i,8); % result : std_u/u*
    result(i,39) = result(i,36)/result(i,8); % result : std_v/u*
    result(i,40) = result(i,37)/result(i,8); % result : std_w/u*
    result(i,41) = result(i,35)/result(i,2); % result : std_u/U
    result(i,42) = result(i,36)/result(i,2); % result : std_v/U
    result(i,43) = result(i,37)/result(i,2); % result : std_w/U
       
    % mean Ts
    result(i,12) = mean(temp(:,4)); % result : mean_Ts 
    Ts = result(i,12); % degree C
%    result(i,48) = skewness(temp(:,4)); % result : skewness of Ts
%    result(i,52) = kurtosis(temp(:,4)); % result : kurtosis of Ts

    % mean ambient T
    result(i,14)=mean(temp(:,8)); % degree C
    amb_tmp = result(i,14); % degree C

    % mean ambient P
    result(i,15)=mean(temp(:,7)); % kPa
    amb_prs = result(i,15); % kPa

    % molar volume of air
    V = 0.082*(273.15+amb_tmp)*10^(-3)/(amb_prs/101.325); % m3/mol   
    
    % mean co2
    result(i,16) = mean(temp(:,5))*V/Mco2; % result : mean co2 
    Xc = result(i,16); % umol/mol (dry air mole fraction)
    result(i,17) = min(temp(:,5))*V/Mco2; % result : min co2
    result(i,18) = max(temp(:,5))*V/Mco2; % result : max co2

    % mean h2o
    result(i,19) = mean(temp(:,6))*V/Mv; % result : mean h2o
    Xv = result(i,19); % mmol/mol (dry air mole fraction)
    result(i,20) = min(temp(:,6))*V/Mv; % result : min h2o
    result(i,21) = max(temp(:,6))*V/Mv; % result : max h2o
    Tk = (Ts + 273.15)/(1+0.000321*Xv); % K ; AmeriFlux
    result(i,13) = Tk-273.15; % result : actual temperature in degree C    
    

    % vapor pressure
    e = Xv*amb_prs/(1000+Xv); % kPa
    result(i,22) = e; % result : vapor pressure in kPa
    es = 6.1365*exp(17.502*(Tk-273.15)/(240.97+(Tk-273.15)))/10; % kPa    
    result(i,23) = es; % result : saturation vapor pressure in kPa
    RH = e/es*100; 
    result(i,24) = RH; % result : relative humidity in percent
    
    % density
    rho_v = Xv*amb_prs*Mv/(R*Tk*(1000+Xv)); % g/m3
    result(i,26) = rho_v; % result : rho_v
%    rho_c = Xc*Mc*10^(-6)*(amb_prs/(R*Tk) - rho_v/Mv); % mgC/m3
    rho_co2 = Xc*10^(-6)*Mco2*(amb_prs/(R*Tk) - rho_v/Mv)*1000; % mg/m3
    result(i,27) = rho_co2; % result : rho_co2
    rho_d = (amb_prs-e)*Md/(R*Tk); % g/m3
    result(i,25) = rho_d; % result : rho_d

    % specific heat capacity
    Cpd = 1005+((Tk-250.03)^2/3364); % J/K/kg ; AmeriFlux
    Cpm = 1859 + 0.13*RH + (Tk-273.15)*(0.193+0.00569*RH) + (Tk-273.15)^2*(0.001+0.000005*RH); % J/K/kg ; AmeriFlux
    rhoCp = Cpm*(rho_v/1000) + Cpd*(rho_d/1000); % specific heat capacity (J/K/m3)
    result(i,28) = rhoCp; % result : rho_Cp (J/K/m3)
    
    % sensible heat flux ; AmeriFlux
    w_bar = result(i,5); % m/s
    q = mean(temp(:,6)); % mean h2o (g/m3)
    cov_wTs = 0;
    cov_wq = 0;
    for j = 1:num_30min_n(i,2)
        cov_wTs = cov_wTs + (temp(j,3)-w_bar)*(temp(j,4)-Ts);
        cov_wq = cov_wq + (temp(j,3)-w_bar)*(temp(j,6)-q); % g/m3 * m/s
    end
    clear j
    cov_wTs = cov_wTs/num_30min_n(i,2);
    cov_wq = cov_wq/num_30min_n(i,2); % g/m2/s
    cov_wq_mol = cov_wq*V/Mv; % mmol*m/s/mol
    Qh = rhoCp * (cov_wTs - (0.000321*Tk*cov_wq_mol)); % AmeriFlux
    result(i,29) = cov_wTs;
    result(i,30) = rhoCp*cov_wTs;
    result(i,31) = Qh;
    
    % latent heat flux ; AmeriFlux
    lambda = 2500.8 - 2.3668*Ts; % J/g
    Qe = lambda*cov_wq; % J/g * g/m2/s => W/m2
    result(i,32) = lambda;
    result(i,33) = Qe; 
    
    % WPL corection for latent heat flux
    u = 1.6077; % ratio of molar masses of air to water
    E = (1+u*(rho_v/rho_d))*(cov_wq + (Qh/rhoCp) * (rho_v/Tk)); % g/m2/s
    Qe_WPL = lambda*E; % W/m2
    result(i,53) = Qe_WPL;
    
    % co2 flux ; AmeriFlux
    co2 = mean(temp(:,5)); % mean co2 (mg/m3)
    cov_wc = 0;
    for j = 1:num_30min_n(i,2)
        cov_wc = cov_wc + (temp(j,3)-w_bar)*(temp(j,5)-co2); % mg/m3 *m/s
    end
    clear j 
    cov_wc = cov_wc/num_30min_n(i,2); % mg/m2/s
    Fc = cov_wc/Mco2; % umol/m2/s
    result(i,34) = Fc;
    
    % WPL corection for co2 flux
    Fc_WPL = cov_wc + u*(E/rho_d)*(rho_co2/(1+u*(rho_v/rho_d))) + ((Qh/rhoCp)*(rho_co2)/Tk); % mg/m2/s
    Fc_WPL = Fc_WPL/Mco2; % umol/m2/s
    result(i,54) = Fc_WPL;

    clear Xv Xc amb_prs
    clear Ts q Tk amb_tmp e Pd rho_v rho_d Cpd es RH Cpm rhoCp w_bar cov_wTs cov_wq_mol cov_wq Qh
    clear lambda V Qe 
    clear co2 cov_wc Fc
    clear E Fc_WPL Qe_WPL rho_co2 u 
end
clear i temp
clear R Mco2 Md Mv
clear po_u po_v po_w po_Ts po_CO2 po_H2O po_amb_prs po_amb_tmp 
clear size_n size_var num_30min_n


total_result{n,1}= dataName(n).name;
total_result{n,2}= result;
clear data result
end

clear sonic_ang n dataName dataDir

toc