%% History
% 130714 Je-Woo Hong, Keunmin Lee and Prof.Jinkyu Hong
% 131001 Je-Woo Hong ; for molar volume, using observed pressure/
%                      for open-path data process (using density)


%% using data
%sonic_ang = 0;
%sonic_ang = 230+8.04; %for SF
%sonic_ang = 220+8.06; %for EP NewTown
sonic_ang = 90+7.24; %for BS Tower

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
        %if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==0.05)        
        if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==0.1)
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

result = zeros(max(num_30min_n(:,1)),34);
% 1; mean wind speed (m/s)
% 2; wind direction (including sonic_angle) (degree)
% 3; u* (m/s)
% 4; std_w/u* (unitless)
% 5; mean Ts (degree C)
% 6; mean CO2 (mg/m3)
% 7; min CO2 (mg/m3)
% 8; max CO2 (mg/m3)
% 9; mean H2O (g/m3)
% 10; min H2O (g/m3)
% 11; max H2O (g/m3)
% 12; mean cell temperature (degree C)
% 13; mean cell pressure (hPa)
% 14; rhoCp (J/K/m3)
% 15; cov_wTs (K*m/s)
% 16; rhoCp*cov_wTs (W/m2)
% 17; H ; AmeriFlux (W/m2) correct the effect of h2o on Ts
% 18; lambda ; latent heat of vaporization of water (J/g)
% 19; LE (W/m2)
% 20; Fco2 (umol/m2/s)
% 21; RH (%)
% 22; e 
% 23; es
% 24; rho_d
% 25; rho_v
% 26; LE_WPL
% 27; Fco2_WPL 
% 28; mean_u (m/s)
% 29; std_v (m/s)
% 30; std_wd (degree)
% 31; n_data (unitless)
% 32; v_bar (m/s)
% 33; w_bar (m/s)
% 34; Tk (degree C)


%% position

po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8;
po_CO2 = 10;
po_H2O = 11;
po_cell_tmpr = 13;
po_cell_prs = 14;


%% 

for i = 1:max(num_30min_n(:,1))
    % extract
    temp = zeros(num_30min_n(i,2),8);
    for j = 1:num_30min_n(i,2)
        if i ==1
            temp(j,1) = data(j,po_u);
            temp(j,2) = data(j,po_v);
            temp(j,3) = data(j,po_w);
            temp(j,4) = data(j,po_Ts);
            temp(j,5) = data(j,po_CO2);
            temp(j,6) = data(j,po_H2O);
            temp(j,7) = data(j,po_cell_tmpr);
            temp(j,8) = data(j,po_cell_prs);
        else
            temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
            temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
            temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
            temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
            temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
            temp(j,6) = data(num_30min_n(i-1,3)+j,po_H2O);
            temp(j,7) = data(num_30min_n(i-1,3)+j,po_cell_tmpr);
            temp(j,8) = data(num_30min_n(i-1,3)+j,po_cell_prs);
        end
    end
    clear j
    
    % mean wind-speed
    temp_ws = 0;
    for j = 1:num_30min_n(i,2)
        temp_ws = temp_ws + (temp(j,1)^2 + temp(j,2)^2 + temp(j,3)^2)^(0.5);
    end
    temp_ws = temp_ws/num_30min_n(i,2);
    result(i,1) = temp_ws;
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
    result(i,30) = temp_std_wd;
    if mean_wd+sonic_ang>360
        mean_wd = mean_wd + sonic_ang-360;
    else
        mean_wd = mean_wd + sonic_ang;
    end
    result(i,2) = mean_wd;
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
    result(i,29)=std(temp(:,2));
    clear u_bar v_bar w_bar alpha beta j temp_wind
    % u* ; u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.25) ; Stull, p.67
    uw = 0; % for mean(u'w')
    vw = 0; % for mean(v'w')
    u_bar = mean(temp(:,1));
    result(i,28) = u_bar;
    v_bar = mean(temp(:,2));
    result(i,32) = v_bar;    
    w_bar = mean(temp(:,3));
    result(i,33) = w_bar;    
    for j = 1:num_30min_n(i,2)
        uw = uw + (temp(j,1)-u_bar)*(temp(j,3)-w_bar);
        vw = vw + (temp(j,2)-v_bar)*(temp(j,3)-w_bar);
    end
    uw = uw/num_30min_n(i,2);
    vw = vw/num_30min_n(i,2);
    result(i,3) = (uw^2 + vw^2)^(0.25);
    clear uw vw u_bar v_bar w_bar j
    
    % std_w/u*
    sigma_w = std(temp(:,3));
    result(i,4) = result(i,3)/sigma_w;
    clear sigma_w;
    
    % mean Ts
    result(i,5) = mean(temp(:,4));
    
    % mean co2
    result(i,6) = mean(temp(:,5));
    result(i,7) = min(temp(:,5));
    result(i,8) = max(temp(:,5));
    
    % mean h2o
    result(i,9) = mean(temp(:,6));
    result(i,10) = min(temp(:,6));
    result(i,11) = max(temp(:,6));
    
    % mean cell T
    result(i,12) = mean(temp(:,7)); % degree C
    
    % mean cell P
    result(i,13) = mean(temp(:,8))*10; % hPa

    % sensible heat flux ; AmeriFlux standard
    Mw = 18.015; % g/mol
    Mc = 44.01; % g/mol
    Ts = mean(temp(:,4)); % degree C
    cell_prs = result(i,13); % hPa    
    cell_tmpr = mean(temp(:,7)); % degree C     
    V = 0.082*(cell_tmpr+273.15)*10^(-3)/(cell_prs/1013.25); % m3/mol  
    q = mean(temp(:,6)); % g/m3
    q_m = mean(temp(:,6))*1000*V/Mw; % mmol/mol
    Tk = (Ts + 273.15)/(1+0.000321*q_m); % K
    e = (q_m/1000)*cell_prs; % hPa
    Pd = cell_prs - e; % hPa
    Rv = 461.495; % J/kg/K
    Rd = 287.05; % J/kg/K
    rho_v = (e*100)/(Rv*(cell_tmpr+273.15)); % kg/m3
    rho_d = (Pd*100)/(Rd*(cell_tmpr+273.15)); % kg/m3
    Cpd = 1005+(((cell_tmpr+273.15)-250.03)^2/3364); % J/K/kg
    es = 6.1365*exp(17.502*cell_tmpr/(240.97+cell_tmpr)); % hPa
    RH = e/es*100; % percent
    Cpm = 1859 + 0.13*RH + (cell_tmpr)*(0.193+0.00569*RH) + (cell_tmpr)^2*(0.001+0.000005*RH); % J/K/kg
    rhoCp = Cpm*rho_v + Cpd*rho_d; % specific heat capacity (J/K/kg)
    w_bar = mean(temp(:,3)); % m/s
    cov_wTs = 0;
    cov_wq = 0;
    for j = 1:num_30min_n(i,2)
        cov_wTs = cov_wTs + (temp(j,3)-w_bar)*(temp(j,4)-Ts);
        cov_wq = cov_wq + (temp(j,3)-w_bar)*(temp(j,6)-q);
    end
    clear j
    cov_wTs = cov_wTs/num_30min_n(i,2);
    cov_wq = cov_wq/num_30min_n(i,2);
    H = rhoCp * (cov_wTs - (0.000321*Tk*cov_wq));
    result(i,14) = rhoCp;
    result(i,15) = cov_wTs;
    result(i,16) = rhoCp*cov_wTs;
    result(i,17) = H;
    
    % latent heat flux ; AmeriFlux Standard
    lambda = 2500.8 - 2.3668*Ts; % J/g
    LE = lambda*cov_wq;
    result(i,18) = lambda;
    result(i,19) = LE;
    
    % co2 flux ; AmeriFlux Standard
    co2 = mean(temp(:,5)); % mg/m3
    cov_wc = 0;
    for j = 1:num_30min_n(i,2)
        cov_wc = cov_wc + (temp(j,3)-w_bar)*(temp(j,5)-co2);
    end
    clear j 
    cov_wc = cov_wc/num_30min_n(i,2); % mg/m2/s
    Fco2 = cov_wc*1000/Mc; % umol/m2/s
    result(i,20) = Fco2; % umol/m2/s
    
    % WPL corection for latent heat flux
    E0 = LE/lambda;
    u = 1.6077;
    E = (1+u*(rho_v/rho_d))*(E0 + (H/rhoCp) * (rho_v/Tk) * 1000); % g/m2/s
    LE_WPL = lambda*E; % W/m2
    
    % WPL corection for co2 flux
    pc = (co2*V/Mc)/1000*cell_prs; % hPa
    Rc = 188.9241; % L*atm/K/mol
    rho_c = (pc*100)/(Rc*Tk); % kg/m3
    Fco2_g = Fco2*Mc/(10^6); % from umolCO2 to gCO2
    Fco2_WPL = Fco2_g+ u*E*(rho_c/rho_d) + (1+u*(rho_v/rho_d))*((1000*rho_c)/Tk)*(H/rhoCp); % gCO2/m2/s
    Fco2_WPL = Fco2_WPL *(10^6)/Mc; % umolCO2/m2/s
    
    result(i,21) = RH;
    result(i,22) = e;
    result(i,23) = es;
    result(i,24) = rho_d;
    result(i,25) = rho_v;
    result(i,26) = LE_WPL;
    result(i,27) = Fco2_WPL;
    result(i,31) = num_30min_n(i,2);
    result(i,34) = Tk - 273.15;
    
    clear Ts q q_m Tk cell_prs e Pd Rv Rd rho_v rho_d Cpd cell_tmpr es RH Cpm rhoCp w_bar cov_wTs cov_wq H 
    clear lambda Mw V LE Mc
    clear co2 cov_wc Fco2
    clear E0 u E LE_WPL
    clear pc Rc rho_c Fco2_g Fco2_WPL
end
clear i temp
clear po_u po_v po_w po_cell_tmpr po_cell_prs po_Ts po_H2O po_CO2 size_n size_var

