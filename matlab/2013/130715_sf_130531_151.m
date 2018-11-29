
data_raw = data;

[size_n, size_var] = size(data);
data = zeros(size_n, size_var+2);

% %doy=151
% start_min = 1440;
% start_sec = 54.8;

% other
start_min = 0000;
start_sec = 0.1;

for i = 1:size_n
    for j = 1:size_var-1
        data(i,j+3) = data_raw(i,j+1);
    end
    if i ==1
        data(i,1) = 151;
        data(i,2) = start_min;
        data(i,3) = start_sec;
    else
        data(i,1) = 151;
        data(i,2) = data(i-1,2);
        data(i,3) = data(i-1,3)+0.1;
        if data(i,3) > 59.95
            data(i,3) = 0.0;
            data(i,2) = data(i,2) +1;
            if mod(data(i,2),100)>59.5
                data(i,2) = data(i,2)+40;
            end
        end
    end
end
clear i j

num_30min = 0;
num_30min_n = zeros(1,3);

for i = 1:size_n
    if i==1
        num_30min = 1;
        num_30min_n(num_30min,1) = num_30min;
        num_30min_n(num_30min,2) = 1;
        num_30min_n(num_30min,3) = 1;
    else
        if ((mod(data(i,2),100)==0)||(mod(data(i,2),100)==30))&&(mod(data(i,3),1)==0.1)
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
% 1; mean wind speed
% 2; wind direction (including sonic_angle)
% 3; u*
% 4; std_w/u*
% 5; mean Ts
% 6; mean CO2 (umol/mol)
% 7; min CO2 (umol/mol)
% 8; max CO2 (umol/mol)
% 9; mean H2O (mmol/mol)
% 10; min H2O (mmol/mol)
% 11; max H2O (mmol/mol)
% 12; mean cell temperature
% 13; mean cell pressure
% 14; rhoCp
% 15; cov_wTs
% 16; rhoCp*cov_wTs
% 17; H ; AmeriFlux
% 18; lambda
% 19; LE
% 20; Fco2
% 21; RH
% 22; e
% 23; es
% 24; rho_d
% 25; rho_v
% 26; LE_WPL
% 27; Fco2_WPL
% 28; mean_u;
% 29; std_v
% 30; std_wd
% 31; n_data
% 32; v_bar
% 33; w_bar
% 34; Tk


%% position

po_u = 4;
po_v = 5;
po_w = 6;
po_Ts = 7;
po_CO2 = 9;
po_H2O = 10;
po_cell_tmpr = 12;
po_cell_prs = 13;


%% 
sonic_ang = 230+8.04; %for SF


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
        mean_wd = mean_wd+sonic_ang-360;
    else
        mean_wd = mean_wd+sonic_ang;
    end
    result(i,2) = mean_wd;
    clear j temp_wd temp_std_wd mean_wd 

    % double rotation
    temp_wind = zeros(num_30min_n(i,2),3);
    alpha = atan(v_bar/u_bar);
    for j = 1:num_30min_n(i,2)
        temp_wind(j,1) = temp(j,1)*cos(alpha)+temp(j,2)*sin(alpha);
        temp_wind(j,2) = -sin(alpha)*temp(j,1)+cos(alpha)*temp(j,2);
        temp_wind(j,3) = temp(j,3);
        temp(j,1) = temp_wind(j,1);
        temp(j,2) = temp_wind(j,2);
        temp(j,3) = temp_wind(j,3);
    end
    u_bar = mean(temp(:,1));
    w_bar = mean(temp(:,3));
    beta = atan(w_bar/u_bar);
    for  j = 1:num_30min_n(i,2)
        temp_wind(j,1) = temp(j,1)*cos(beta)+temp(j,3)*sin(beta);
        temp_wind(j,2) = temp(j,2);
        temp_wind(j,3) = -1*temp(j,1)*sin(beta) + temp(j,3)*cos(beta);
        temp(j,1) = temp_wind(j,1);
        temp(j,2) = temp_wind(j,2);
        temp(j,3) = temp_wind(j,3);
    end
    result(i,29)=std(temp(:,2));
    result(i,32)=v_bar;
    result(i,33)=w_bar;
    clear u_bar v_bar w_bar alpha beta j temp_wind

    % u* ; u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.5) ; Stull, p.67
    uw = 0; % for mean(u'w')
    vw = 0; % for mean(v'w')
    u_bar = mean(temp(:,1));
    result(i,28) = u_bar;
    v_bar = mean(temp(:,2));
    w_bar = mean(temp(:,3));
    for j = 1:num_30min_n(i,2)
        uw = uw + (temp(j,1)-u_bar)*(temp(j,3)-w_bar);
        vw = vw + (temp(j,2)-v_bar)*(temp(j,3)-w_bar);
    end
    uw = uw/num_30min_n(i,2);
    vw = vw/num_30min_n(i,2);
    result(i,3) = (uw^2 + vw^2)^(0.5);
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
    result(i,12) = mean(temp(:,7));
    
    % mean cell P
    result(i,13) = mean(temp(:,8));

    % sensible heat flux ; AmeriFlux standard
    Ts = mean(temp(:,4));
    q = mean(temp(:,6));
    Tk = (Ts + 273.15)/(1+0.000321*q);
    cell_prs = mean(temp(:,8));
    e = (q*10^(-3))*cell_prs;
    Pd = cell_prs - e;
    Rv = 461.495;
    Rd = 287.05;
    cell_tmpr = mean(temp(:,7));
    rho_v = (e*1000)/(Rv*(cell_tmpr+273.15));
    rho_d = (Pd*1000)/(Rd*(cell_tmpr+273.15));
    Cpd = 1005+(((cell_tmpr+273.15)-250.03)^2/3364);
    es = 6.1365*exp(17.502*cell_tmpr/(240.97+cell_tmpr))*0.1;
    RH = e/es*100;
    Cpm = 1859 + 0.13*RH + (cell_tmpr)*(0.193+0.00569*RH) + (cell_tmpr)^2*(0.001+0.000005*RH);
    rhoCp = Cpm*rho_v + Cpd*rho_d;
    w_bar = mean(temp(:,3));
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
    lambda = 2500.8 - 2.3668*Ts;
    Mw = 18.015;
    V = 0.082*Tk/1.10325;
    LE = lambda*Mw*cov_wq/(V);
    result(i,18) = lambda;
    result(i,19) = LE;
    
    % co2 flux ; AmeriFlux Standard
    co2 = mean(temp(:,5));
    cov_wc = 0;
    for j = 1:num_30min_n(i,2)
        cov_wc = cov_wc + (temp(j,3)-w_bar)*(temp(j,5)-co2);
    end
    clear j 
    cov_wc = cov_wc/num_30min_n(i,2);
    Fco2 = cov_wc/(10^(-3)*V);
    result(i,20) = Fco2;
    
    % WPL corection for latent heat flux
    E0 = LE/lambda;
    u = 1.6077;
    E = (1+u*(rho_v/rho_d))*(E0 + (H/rhoCp) * (rho_v/Tk));
    LE_WPL = lambda*E;
    
    % WPL corection for co2 flux
    pc = (co2*10^(-6))*cell_prs;
    Rc = 188.9241; % from  
    rho_c = (pc*1000)/(Rc*(cell_tmpr+273.15));
    Fco2_kg = Fco2 * 12/(10^6);
    Fco2_WPL = Fco2_kg + u*(E/rho_d)*(rho_c/(1+u*(rho_v/rho_d)))+(H/rhoCp)*(rho_c/Tk);
    Fco2_WPL = Fco2_WPL *(10^6)/12;
    
    result(i,21) = RH;
    result(i,22) = e;
    result(i,23) = es;
    result(i,24) = rho_d;
    result(i,25) = rho_v;
    result(i,26) = LE_WPL;
    result(i,27) = Fco2_WPL;
    result(i,31) = num_30min_n(i,2);
    
    clear Ts q Tk cell_prs e Pd Rv Rd rho_v rho_d Cpd cell_tmpr es RH Cpm rhoCp w_bar cov_wTs cov_wq H
    clear lambda Mw V LE 
    clear co2 cov_wc Fco2
    clear E0 u E LE_WPL
    clear pc Rc rho_c Fco2_kg Fco2_WPL
end
clear i temp
clear po_u po_v po_w po_cell_tmpr po_cell_prs po_Ts po_H2O po_CO2 size_n size_var

