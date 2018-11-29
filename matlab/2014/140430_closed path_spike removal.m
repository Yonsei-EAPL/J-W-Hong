%% History
% 130714 Mr.Je-Woo Hong & Prof.Jinkyu Hong, discussed with Mr.Keunmin Lee 
% 131023 Mr.Je-Woo Hong; Removal of WPL terms (closed-path)
% 140430 Mr.Je-Woo Hong & Prof.Jinkyu Hong, discussed with Mr.K-H Hong (B&P) 
    % about amb_press and amb_temperature


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
% 14; mean cell pressure (kPa)
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
% 27; rho_c (mgC/m3)
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
    result(i,46)=skewness(temp(:,1)); % result : skewness of u
    result(i,50)=kurtosis(temp(:,1)); % result : kurtosis of u
    result(i,47)=skewness(temp(:,1)); % result : skewness of u
    result(i,51)=kurtosis(temp(:,1)); % result : kurtosis of u    
    clear u_bar v_bar w_bar alpha beta j temp_wind

    % friction velocity
        % u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.5) 
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
    result(i,8) = (uw^2 + vw^2)^(0.5); % result : u_star
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
    result(i,48) = skewness(temp(:,4)); % result : skewness of Ts
    result(i,52) = kurtosis(temp(:,4)); % result : kurtosis of Ts
        
    % mean co2
    result(i,16) = mean(temp(:,5)); % result : mean co2 
    Xc = result(i,16); % umol/mol (dry air mole fraction)
    result(i,17) = min(temp(:,5)); % result : min co2
    result(i,18) = max(temp(:,5)); % result : max co2
    
    % mean h2o
    result(i,19) = mean(temp(:,6)); % result : mean h2o
    Xv = result(i,19); % mmol/mol (dry air mole fraction)
    result(i,20) = min(temp(:,6)); % result : min h2o
    result(i,21) = max(temp(:,6)); % result : max h2o
    Tk = (Ts + 273.15)/(1+0.000321*Xv); % K ; AmeriFlux
    result(i,13) = Tk-273.15; % result : actual temperature in degree C    
    
    % mean cell P
    result(i,14) = mean(temp(:,7)); % kPa
    
    % mean ambient P
    result(i,15) = (mean(temp(:,7))-mean(temp(:,8))); % kPa    
    amb_prs = result(i,15); % kPa
    
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
    rho_c = Xc*Mc*10^(-6)*(amb_prs/(R*Tk) - rho_v/Mv); % mgC/m3
    result(i,27) = rho_c; % result : rho_c
    rho_d = (amb_prs-e)*Md/(R*Tk); % g/m3
    result(i,25) = rho_d; % result : rho_d

    % specific heat capacity
    Cpd = 1005+((Tk-250.03)^2/3364); % J/K/kg ; AmeriFlux
    Cpm = 1859 + 0.13*RH + (Tk-273.15)*(0.193+0.00569*RH) + (Tk-273.15)^2*(0.001+0.000005*RH); % J/K/kg ; AmeriFlux
    rhoCp = Cpm*(rho_v/1000) + Cpd*(rho_d/1000); % specific heat capacity (J/K/kg)
    result(i,28) = rhoCp; % result : rho_Cp (J/K/kg)
    
    % sensible heat flux ; AmeriFlux
    w_bar = result(i,5); % m/s
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
    result(i,29) = cov_wTs;
    result(i,30) = rhoCp*cov_wTs;
    result(i,31) = Qh;
    
    % latent heat flux ; AmeriFlux
    lambda = 2500.8 - 2.3668*Ts; % J/g
    V = 0.082*Tk*10^(-3)/(amb_prs/101.325); % m3/mol
    Qe = lambda*Mv*cov_wq/V; 
    result(i,32) = lambda;
    result(i,33) = Qe; 
    
    % co2 flux ; AmeriFlux
    co2 = result(i,16); % mean co2 (umol/mol)
    cov_wc = 0;
    for j = 1:num_30min_n(i,2)
        cov_wc = cov_wc + (temp(j,3)-w_bar)*(temp(j,5)-co2);
    end
    clear j 
    cov_wc = cov_wc/num_30min_n(i,2); % umol/mol * m/s
    Fc = cov_wc/V; % umol/m2/s 
    result(i,34) = Fc;

    clear Xv Xc amb_prs
    clear Ts q Tk cell_prs e Pd rho_v rho_d rho_c Cpd es RH Cpm rhoCp w_bar cov_wTs cov_wq Qh
    clear lambda V Qe 
    clear co2 cov_wc Fc
end
clear i temp
clear R Mc Md Mv
clear po_u po_v po_w po_Ts po_CO2 po_H2O po_cell_prs po_diff_prs 
clear size_n size_var sonic_ang num_30min_n


