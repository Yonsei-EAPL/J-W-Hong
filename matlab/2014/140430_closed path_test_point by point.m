%% History
% 130714 Mr.Je-Woo Hong & Prof.Jinkyu Hong, discussed with Mr.Keunmin Lee 
% 131023 Mr.Je-Woo Hong; Removal of WPL terms (closed-path)
% 140430 Mr.Je-Woo Hong & Prof.Jinkyu Hong, discussed with Mr.K-H Hong (B&P) 
    % about amb_press and amb_temperature
% 140502 Mr.Je-Woo Hong & Prof.Jinkyu Hong
    % for test WPL terms for closed-path


%% input(*.csv) information

% sonic_ang = 230+8.04; %for SF
sonic_ang = 220+8.06; %for EP NewTown

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

result = zeros(max(num_30min_n(:,1)),52);
% 1 ; n_data (unitless)
% 2 ; mean Ts (degree C)
% 3 ; mean Tk (degree C)
% 4 ; mean cell pressure (kPa)
% 5 ; mean ambient pressure (kPa)
% 6 ; mean CO2 (umol/mol)
% 7 ; mean H2O (mmol/mol)
% 8 ; Qh ; AmeriFlux (W/m2) correct the effect of h2o on Ts
% 9 ; Qe (W/m2)
% 10; Fc (umol/m2/s)
% 11; Fc2 (umol/m2/s)


%% position

po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8;
po_CO2 = 10;
po_H2O = 11;
po_cell_prs = 14;
po_diff_prs = 17;


%% constant
R = 8.3143*10^(-6); % kPa*m3/K/mmol
Mc = 12; % mgC/mmol
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
            temp(j,7) = data(j,po_cell_prs);
            temp(j,8) = data(j,po_diff_prs);
        else
            temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
            temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
            temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
            temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
            temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
            temp(j,6) = data(num_30min_n(i-1,3)+j,po_H2O);
            temp(j,7) = data(num_30min_n(i-1,3)+j,po_cell_prs);
            temp(j,8) = data(num_30min_n(i-1,3)+j,po_diff_prs);
        end
    end
    clear j
    
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
    clear u_bar v_bar w_bar alpha beta j temp_wind

    % mean Ts
    result(i,2) = mean(temp(:,4)); % result : mean_Ts 
    Ts = result(i,2); % degree C
    % mean co2
    result(i,6) = mean(temp(:,5)); % result : mean co2 
    Xc = result(i,6); % umol/mol (dry air mole fraction)
    % mean h2o
    result(i,7) = mean(temp(:,6)); % result : mean h2o
    Xv = result(i,7); % mmol/mol (dry air mole fraction)
    Tk = (Ts + 273.15)/(1+0.000321*Xv); % K ; AmeriFlux
    result(i,3) = Tk-273.15; % result : actual temperature in degree C    
    % mean cell P
    result(i,4) = mean(temp(:,7)); % kPa
    % mean ambient P
    result(i,5) = (mean(temp(:,7))-mean(temp(:,8))); % kPa    
    amb_prs = result(i,5); % kPa
    % vapor pressure
    e = Xv*amb_prs/(1000+Xv); % kPa
    es = 6.1365*exp(17.502*(Tk-273.15)/(240.97+(Tk-273.15)))/10; % kPa    
    RH = e/es*100; 
    % density
    rho_v = Xv*amb_prs*Mv/(R*Tk*(1000+Xv)); % g/m3
    rho_c = Xc*Mc*10^(-6)*(amb_prs/(R*Tk) - rho_v/Mv); % mgC/m3
    rho_d = (amb_prs-e)*Md/(R*Tk); % g/m3
    % specific heat capacity
    Cpd = 1005+((Tk-250.03)^2/3364); % J/K/kg ; AmeriFlux
    Cpm = 1859 + 0.13*RH + (Tk-273.15)*(0.193+0.00569*RH) + (Tk-273.15)^2*(0.001+0.000005*RH); % J/K/kg ; AmeriFlux
    rhoCp = Cpm*(rho_v/1000) + Cpd*(rho_d/1000); % specific heat capacity (J/K/kg)
    % sensible heat flux ; AmeriFlux
    w_bar = mean(temp(:,3)); % m/s
    cov_wTs = 0;
    cov_wq = 0;
    for j = 1:num_30min_n(i,2)
        cov_wTs = cov_wTs + (temp(j,3)-w_bar)*(temp(j,4)-Ts);
        cov_wq = cov_wq + (temp(j,3)-w_bar)*(temp(j,6)-Xv); % mmol/mol *m/s
    end
    clear j
    cov_wTs = cov_wTs/num_30min_n(i,2);
    cov_wq = cov_wq/num_30min_n(i,2);
    Qh = rhoCp * (cov_wTs - (0.000321*Tk*cov_wq)); % AmeriFlux
    result(i,8) = Qh;
    % latent heat flux ; AmeriFlux
    lambda = 2500.8 - 2.3668*Ts; % J/g
    V = 0.082*Tk*10^(-3)/(amb_prs/101.325); % m3/mol
    Qe = lambda*Mv*cov_wq/V; 
    result(i,9) = Qe; 
    % co2 flux ; AmeriFlux
    co2 = result(i,6); % mean co2 (umol/mol)
    cov_wc = 0;
    for j = 1:num_30min_n(i,2)
        cov_wc = cov_wc + (temp(j,3)-w_bar)*(temp(j,5)-co2);
    end
    clear j 
    cov_wc = cov_wc/num_30min_n(i,2); % umol/mol * m/s
    Fc = cov_wc/V; % umol/m2/s 
    result(i,10) = Fc;
    
    
    
    clear Xv Xc amb_prs
    clear Ts q Tk cell_prs e Pd rho_v rho_d rho_c Cpd es RH Cpm rhoCp w_bar cov_wTs cov_wq Qh
    clear lambda V Qe 
    clear co2 cov_wc Fc
end
clear i temp
clear R Mc Md Mv
clear po_u po_v po_w po_Ts po_CO2 po_H2O po_cell_prs po_diff_prs 
clear size_n size_var sonic_ang num_30min_n


