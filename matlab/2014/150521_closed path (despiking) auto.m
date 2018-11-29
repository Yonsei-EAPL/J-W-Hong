%% History
% 130714 Mr.Je-Woo Hong & Prof.Jinkyu Hong, discussed with Mr.Keunmin Lee 
% 131023 Mr.Je-Woo Hong
    % remove the WPL terms (closed-path)
% 140430 Mr.Je-Woo Hong & Prof.Jinkyu Hong, discussed with Mr.K-H Hong (B&P) 
    % about amb_press and amb_temperature
% 150115 Mr.Keunmin Lee
    % kurtosis, and skewness modify (from "u, u, u" to "u, v, w")
% 150521 Mr.Je-Woo Hong & Keunmin Lee
    % TKE, and despiking process
    % automatic run in specific directory
% 150605 Mr.Je-Woo Hong
    % diagnostic value, stationary test    
    

%% position

po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8;
po_diag1 = 9; %% 2015-06-05
po_CO2 = 10;
po_H2O = 11;
po_diag2 = 12; %% 2015-06-05
po_cell_prs = 14;
po_diff_prs = 17;


%% constant
R = 8.3143*10^(-6); % kPa*m3/K/mmol
Mc = 12; % mgC/mmol
Md = 0.029; % g/mmol
Mv = 0.018015; % g/mmol
des_sigma = 3.5; % thresohold standard deviation value for despiking
des_window = 5; % window width in minute
    des_window = des_window * 60 * 10;
% sonic_ang = 230+8.04; %for SF
sonic_ang = 220+8.06; %for EP NewTown


%% input(*.csv) information
dataDir = 'E:\EAPL\JW_Observation\은평뉴타운\_csv'; % folder name
dataName = dir(fullfile(dataDir, 'CSV_3302.ts_data_*'))'; 
total_result=cell(length(dataName),2);

wait = waitbar(0,'...Flux Data Process...');

n_st = 141;
n_en = 150;
for t=n_st:n_en%1:length(dataName) % 여러개의 CSV 파일을 연속적으로 처리하기 위한 loop

    data = importdata(fullfile(dataDir, dataName(t).name));  

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

    result = zeros(max(num_30min_n(:,1)),62);
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
    % 53; TKE (m2/s2)
    % 54; n_data for diagnostic value(CSAT3)
    % 55; n_data for diagnostic value(EC155)
    % 56; stationary test for w'Ts'
    % 57; stationary test for w'q'
    % 58; stationary test for Fc
    % 59; n_data after despiking 
    % 60; Qh after despiking
    % 61; Qe after despiking
    % 62; Fc after despiking


    %% main process
    for i = 1:max(num_30min_n(:,1))
        waitbar(((t-n_st)/(n_en-n_st)),wait,sprintf('%3f',i))
        %waitbar((t/length(dataName)),wait,sprintf('%3f',i))
        result(i,1) = num_30min_n(i,2);    
        if num_30min_n(i,2)==18000
            % extract
            temp = zeros(num_30min_n(i,2),10);
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
                    temp(j,9) = data(j,po_diag1);
                    temp(j,10) = data(j,po_diag2);                    
                else
                    temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
                    temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
                    temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
                    temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
                    temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
                    temp(j,6) = data(num_30min_n(i-1,3)+j,po_H2O);
                    temp(j,7) = data(num_30min_n(i-1,3)+j,po_cell_prs);
                    temp(j,8) = data(num_30min_n(i-1,3)+j,po_diff_prs);
                    temp(j,9) = data(num_30min_n(i-1,3)+j,po_diag1);
                    temp(j,10) = data(num_30min_n(i-1,3)+j,po_diag2);
                end
            end
            clear j

            % check diagnostic value : 2015-06-04
            temp_diag_n1=0;
            temp_diag_n2=0;
            for j = 1:18000
                if temp(j,9)~=0
                    temp_diag_n1 = temp_diag_n1+1;
                end
                if temp(j,10)~=0
                    temp_diag_n2 = temp_diag_n2+1;
                end
            end
            if (temp_diag_n1>0)&&(temp_diag_n1<18000)
                for j = 1:18000
                    if temp(j,9)~=0
                        if j==1
                            if (temp(j+1,9)==0)&&(temp(j+2,9)==0)
                                temp(j,1) = temp(j+1,1)*2-temp(j+2,1);
                                temp(j,2) = temp(j+1,2)*2-temp(j+2,2);
                                temp(j,3) = temp(j+1,3)*2-temp(j+2,3);
                                temp(j,4) = temp(j+1,4)*2-temp(j+2,4);
                            else
                                k=0;
                                while(temp(j+k,9)~=0)
                                    k=k+1;
                                end
                                temp(j,1) = temp(j+k,1);
                                temp(j,2) = temp(j+k,2);
                                temp(j,3) = temp(j+k,3);
                                temp(j,4) = temp(j+k,4);
                                clear k
                            end
                        elseif j==18000
                            if (temp(j-1,9)==0)&&(temp(j-2,9)==0)                    
                                temp(j,1) = temp(j-1,1)*2-temp(j-2,1);
                                temp(j,2) = temp(j-1,2)*2-temp(j-2,2);
                                temp(j,3) = temp(j-1,3)*2-temp(j-2,3);
                                temp(j,4) = temp(j-1,4)*2-temp(j-2,4);
                            else
                                k=0;
                                while(temp(j-k,9)~=0)
                                    k=k+1;
                                end
                                temp(j,1) = temp(j-k,1);
                                temp(j,2) = temp(j-k,2);
                                temp(j,3) = temp(j-k,3);
                                temp(j,4) = temp(j-k,4);
                                clear k                        
                            end
                        else
                            if (temp(j-1,9)==0)&&(temp(j+1,9)==0)                                        
                                temp(j,1) = (temp(j-1,1)+temp(j+1,1))/2;
                                temp(j,2) = (temp(j-1,2)+temp(j+1,2))/2;
                                temp(j,3) = (temp(j-1,3)+temp(j+1,3))/2;
                                temp(j,4) = (temp(j-1,4)+temp(j+1,4))/2;
                            else
                                k=0;
                                while((temp(j+k,9)~=0)&&(j+k<18000))
                                    k=k+1;
                                end
                                temp(j,1) = temp(j+k,1);
                                temp(j,2) = temp(j+k,2);
                                temp(j,3) = temp(j+k,3);
                                temp(j,4) = temp(j+k,4);
                                clear k
                            end
                        end
                    end            
                end
            end
            if (temp_diag_n2>0)&&(temp_diag_n2<18000)
                for j = 1:18000
                    if temp(i,10)~=0
                        if j==1
                            if (temp(j+1,10)==0)&&(temp(j+2,10)==0)
                                temp(j,5) = temp(j+1,5)*2-temp(j+2,5);
                                temp(j,6) = temp(j+1,6)*2-temp(j+2,6);
                            else
                                k=0;
                                while(temp(j+k,10)~=0)
                                    k=k+1;
                                end
                                temp(j,5) = temp(j+k,5);
                                temp(j,6) = temp(j+k,6);
                                clear k
                            end
                        elseif j==18000
                            if (temp(j-1,10)==0)&&(temp(j-2,10)==0)
                                temp(j,5) = temp(j-1,5)*2-temp(j-2,5);
                                temp(j,6) = temp(j-1,6)*2-temp(j-2,6);
                            else
                                k=0;
                                while(temp(j-k,10)~=0)
                                    k=k+1;
                                end
                                temp(j,5) = temp(j-k,5);
                                temp(j,6) = temp(j-k,6);
                                clear k
                            end
                        else
                            if (temp(j-1,10)==0)&&(temp(j+1,10)==0)
                                temp(j,5) = (temp(j-1,5)+temp(j+1,5))/2;
                                temp(j,6) = (temp(j-1,6)+temp(j+1,6))/2;
                            else
                                k=0;
                                while((temp(j+k,10)~=0)&&(j+k<18000))
                                    k=k+1;
                                end
                                temp(j,5) = temp(j+k,5);
                                temp(j,6) = temp(j+k,6);
                                clear k
                            end
                        end
                    end
                end
            end 
            clear j temp_diag
            result(i,54) = temp_diag_n1; % n_data for diagnostic value(CSAT)
            result(i,55) = temp_diag_n2; % n_data for diagnostic value(EC155)
            clear temp_diag_n1 temp_diag_n2
            
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
                % u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.5) 
                % Stull, p.67
            uw = 0; % for mean(u'w')
            vw = 0; % for mean(v'w')
            TKE = 0; % for TKE
            u_bar = mean(temp(:,1));
            result(i,3) = u_bar; % result : mean u
            v_bar = mean(temp(:,2));
            result(i,4) = v_bar; % result : mean v
            w_bar = mean(temp(:,3));
            result(i,5) = w_bar; % result : mean w
            for j = 1:num_30min_n(i,2)
                uw = uw + (temp(j,1)-u_bar)*(temp(j,3)-w_bar);
                vw = vw + (temp(j,2)-v_bar)*(temp(j,3)-w_bar);
                TKE = TKE + (temp(j,1)-u_bar)^2+(temp(j,3)-w_bar)^2+(temp(j,2)-v_bar)^2;
            end
            uw = uw/num_30min_n(i,2);
            vw = vw/num_30min_n(i,2);
            TKE = TKE/2/num_30min_n(i,2);
            result(i,8) = (uw^2 + vw^2)^(0.5); % result : u_star
            result(i,9) = uw; % result : cov_uw
            result(i,44) = uw/result(i,35)/result(i,37);
            result(i,53) = TKE; % result : TKE
            clear uw vw u_bar v_bar w_bar TKE j

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
        % stationary test
        for j = 1:6
            cov_wTs =0;
            cov_wq =0;
            cov_wc =0;
            temp_wind = zeros(3000,3);
            u_bar2 = mean(temp((j-1)*3000+1:j*3000,1));
            v_bar2 = mean(temp((j-1)*3000+1:j*3000,2));
            w_bar2 = mean(temp((j-1)*3000+1:j*3000,3));
            alpha = atan2(v_bar2,u_bar2);
            beta = atan2(w_bar2,((u_bar2^2 + v_bar2^2)^(0.5)));
            for k = 1:3000
                temp_wind(k,1) = cos(beta)*(cos(alpha)*temp((j-1)*3000+k,1)+sin(alpha)*temp((j-1)*3000+k,2))+sin(beta)*temp((j-1)*3000+k,3);
                temp_wind(k,2) = -sin(alpha)*temp((j-1)*3000+k,1)+cos(alpha)*temp((j-1)*3000+k,2);
                temp_wind(k,3) = -sin(beta)*(cos(alpha)*temp((j-1)*3000+k,1)+sin(alpha)*temp((j-1)*3000+k,2))+cos(beta)*temp((j-1)*3000+k,3);
            end
            clear k
            w_bar2 = mean(temp_wind(:,3));
            Ts2 = mean(temp((j-1)*3000+1:j*3000,4));
            Xv2 = mean(temp((j-1)*3000+1:j*3000,6));
            co22 = mean(temp((j-1)*3000+1:j*3000,5));
            for k = 1:3000
                cov_wTs = cov_wTs + (temp_wind(k,3)-w_bar2)*(temp(k,4)-Ts2);
                cov_wq = cov_wq + (temp_wind(k,3)-w_bar2)*(temp(k,6)-Xv2); % mmol/mol *m/s
                cov_wc = cov_wc + (temp_wind(k,3)-w_bar2)*(temp(k,5)-co22);
            end
            clear k
            result(i,56) =result(i,56)+ cov_wTs/3000;
            result(i,57) =result(i,57)+ cov_wq/3000;
            result(i,58) =result(i,58)+ cov_wc/3000;
        end
        clear j
        clear u_bar2 v_bar2 w_bar2 Ts2 Xv2 co22 cov_wTs cov_wq cov_wc temp_wind alpha beta
        result(i,56) =result(i,56)/6;
        result(i,57) =result(i,57)/6;
        result(i,58) =result(i,58)/6/V;
        result(i,56) = rhoCp * (result(i,56) - (0.000321*Tk*result(i,57)));
        result(i,57) =lambda*Mv*result(i,57)/V;
        
        % despiking
        des_n = 0;
        for j = 1:(18000/des_window)*5-4
            des_n_temp = 0;
            des_n_before = des_n;
            for n =1:10
                a = (j-1)*600+1-des_n_before;
                b = (j-1)*600+des_window-des_n_temp-des_n_before;
                des_std = abs(std(temp(a:b,5)));
                des_mean = mean(temp(a:b,5));
                for k = 1:des_window-des_n_temp
                    if ((j-1)*600+k)<(18000-des_n)
                        if abs(temp((j-1)*600-des_n_before+k,5)-des_mean)>(des_sigma+(n-1)*0.1)*des_std
                            temp((j-1)*600-des_n_before+k,:)=[];
                            des_n = des_n +1;
                            des_n_temp = des_n_temp +1;
                        end
                    end
                end
            end
        end
        [a b] = size(temp);
        result(i,59) = a;
        clear a b des_n des_std des_mean 

        % after despiking
        u_bar = mean(temp(:,1));
        v_bar = mean(temp(:,2));
        w_bar = mean(temp(:,3));
        Ts = mean(temp(:,4)); % result : mean_Ts 
        Xc = mean(temp(:,5)); % result : mean co2 
        Xv = mean(temp(:,6)); % result : mean h2o
        Tk = (Ts + 273.15)/(1+0.000321*Xv); % K ; AmeriFlux
        amb_prs = (mean(temp(:,7))-mean(temp(:,8))); % kPa    
        e = Xv*amb_prs/(1000+Xv); % kPa
        es = 6.1365*exp(17.502*(Tk-273.15)/(240.97+(Tk-273.15)))/10; % kPa    
        RH = e/es*100; 
        rho_v = Xv*amb_prs*Mv/(R*Tk*(1000+Xv)); % g/m3
        rho_c = Xc*Mc*10^(-6)*(amb_prs/(R*Tk) - rho_v/Mv); % mgC/m3
        rho_d = (amb_prs-e)*Md/(R*Tk); % g/m3
        Cpd = 1005+((Tk-250.03)^2/3364); % J/K/kg ; AmeriFlux
        Cpm = 1859+0.13*RH+(Tk-273.15)*(0.193+0.00569*RH)+(Tk-273.15)^2*(0.001+0.000005*RH); % J/K/kg ; AmeriFlux
        rhoCp = Cpm*(rho_v/1000)+Cpd*(rho_d/1000); % specific heat capacity (J/K/kg)
        cov_wTs = 0;
        cov_wq = 0;
        for j = 1:result(i,59)
            cov_wTs = cov_wTs + (temp(j,3)-w_bar)*(temp(j,4)-Ts);
            cov_wq = cov_wq + (temp(j,3)-w_bar)*(temp(j,6)-Xv); % mmol/mol *m/s
        end
        clear j
        cov_wTs = cov_wTs/result(i,59);
        cov_wq = cov_wq/result(i,59);
        Qh = rhoCp*(cov_wTs-(0.000321*Tk*cov_wq)); % AmeriFlux
        result(i,60) = Qh;
        lambda = 2500.8-2.3668*Ts; % J/g
        V = 0.082*Tk*10^(-3)/(amb_prs/101.325); % m3/mol
        Qe = lambda*Mv*cov_wq/V; 
        result(i,61) = Qe; 
        co2 = Xc; % mean co2 (umol/mol)
        cov_wc = 0;
        for j = 1:result(i,59)
            cov_wc = cov_wc + (temp(j,3)-w_bar)*(temp(j,5)-co2);
        end
        clear j 
        cov_wc = cov_wc/result(i,59); % umol/mol * m/s
        Fc = cov_wc/V; % umol/m2/s 
        result(i,62) = Fc;

        clear Xv Xc amb_prs
        clear Ts q Tk cell_prs e Pd rho_v rho_d rho_c Cpd es RH Cpm rhoCp w_bar cov_wTs cov_wq Qh
        clear lambda V Qe 
        clear co2 cov_wc Fc
        clear u_bar v_bar 
        end
    end
    total_result{t,1}= dataName(t).name;
    total_result{t,2}= result;
    
    clear i
    clear size_n size_var num_30min_n
end
close(wait);
clear R Mc Md Mv
clear po_u po_v po_w po_Ts po_CO2 po_H2O po_cell_prs po_diff_prs sonic_ang 



