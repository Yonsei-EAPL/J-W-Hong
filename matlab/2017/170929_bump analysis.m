%% History
%2017-Sep-29 JW Hong : for bump analysis


%% information
% site : EP
% IRGA : EC-155
% period :


%% site-period dependent information (location)
% EP-
sonic_ang = 220-8.06; %2016-Sep-21 KM (minus!! not plus)
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


%% input(*.csv) information
dataDir = 'E:\EAPL\JW_Observation\은평뉴타운\_csv'; % folder name % EP
dataName = dir(fullfile(dataDir, 'CSV_3302.ts_data_2016_03*'))'; % EP
total_result=cell(length(dataName),2);


%% main
% n_st = 141;
% n_en = 150;
for t=1:length(dataName) %n_st:n_en %loop for all data-files(*.csv)
    data = importdata(fullfile(dataDir, dataName(t).name));  % import csv file
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
    
    result = zeros(max(num_30min_n(:,1)),62); % define: result
    % 1; n_data (unitless)
    % 2; mean wind speed, U (m/s)
    % 3; u_bar (m/s)
    % 4; v_bar (m/s)
    % 5; w_bar (m/s)
    % 6; wind direction (including sonic_angle) (degree)
    % 7; std_wd (degree)
    % 8; u* (m/s)
    % 9; co2 flux; Fc (umol/m2/s) original
    %10; co2 flux; Fc (umol/m2/s) without bump signal
    %11; co2 flux; Fc (umol/m2/s) only bumpsignal
    %12; var_co2 (umol2/mol2) variance of co2 original
    %13; var_co2 (umol2/mol2) variance of co2 without bump signal
    %14; var_co2 (umol2/mol2) variance of co2 only bumpsignal
    
    
    for i = 1:max(num_30min_n(:,1))
        result(i,1) = num_30min_n(i,2);
        if num_30min_n(i,2)==18000 % if n is less/more than 18000 then skip the process
            temp = zeros(num_30min_n(i,2),10); % extract 30-min raw data
            for j = 1:num_30min_n(i,2)
                if i == 1
                    temp(j,1) = data(j,po_u); % u
                    temp(j,2) = data(j,po_v); % v
                    temp(j,3) = data(j,po_w); % w
                    temp(j,4) = data(j,po_Ts); % Ts
                    temp(j,5) = data(j,po_CO2); % co2
                    temp(j,6) = data(j,po_H2O); % h2o
                    temp(j,7) = data(j,po_cell_prs); % cell_pressure
                    temp(j,8) = data(j,po_diff_prs); % diff
                    %temp(j,9) = data(j,po_diag1); no use diag; co2 wo bump
                    %temp(j,10) = data(j,po_diag2); no use diag; bump only
                else
                    temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
                    temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
                    temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
                    temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
                    temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
                    temp(j,6) = data(num_30min_n(i-1,3)+j,po_H2O);
                    temp(j,7) = data(num_30min_n(i-1,3)+j,po_cell_prs);
                    temp(j,8) = data(num_30min_n(i-1,3)+j,po_diff_prs);
                    %temp(j,9) = data(num_30min_n(i-1,3)+j,po_diag1); no use diag; co2 wo bump
                    %temp(j,10) = data(num_30min_n(i-1,3)+j,po_diag2); no use diag; bump only
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
            
            % mean wind-direction (51 lines)
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
            beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.25)));
            for j = 1:num_30min_n(i,2)
                temp_wind(j,1) = cos(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+sin(beta)*temp(j,3);
                temp_wind(j,2) = -sin(alpha)*temp(j,1)+cos(alpha)*temp(j,2);
                temp_wind(j,3) = -sin(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+cos(beta)*temp(j,3);
                temp(j,1) = temp_wind(j,1); % rotated u
                temp(j,2) = temp_wind(j,2); % rotated v
                temp(j,3) = temp_wind(j,3); % rotated w
            end
            clear u_bar v_bar w_bar alpha beta j temp_wind            
            
            % friction velocity % u*=(mean(u'w')^2+mean(v'w')^2)^(0.25) % Stull, p.67
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
            clear uw vw u_bar v_bar j

            Ts = mean(temp(:,4)); % degreeC
            Xc = mean(temp(:,5)); % umol/mol (dry air mole fraction)
            Xv = mean(temp(:,6)); % mmol/mol (dry air mole fraction)
            Tk = (Ts + 273.15)/(1+0.000321*Xv) % Kelvin ; AmeriFlux
            amb_prs = (mean(temp(:,7))-mean(temp(:,8))); % mean ambient P (=cellP - diffP) % kPa
            V = 0.082*Tk*10^(-3)/(amb_prs/101.325); % m3/mol
            clear Ts Tk Xv amb_prs
            
            cov_wc = 0; % co2 flux (Fc)
            for j = 1:num_30min_n(i,2)
                cov_wc = cov_wc + (temp(j,3)-w_bar)*(temp(j,5)-Xc);
            end
            clear j
            cov_wc = cov_wc/num_30min_n(i,2); % umol/mol * m/s
            Fc = cov_wc/V; % original co2 flux ; umol/m2/s
            result(i,9) = Fc; % result : co2 flux original
            result(i,12) = var(temp(:,5)); % result : var_co2 original
            clear Xc cov_wc Fc
            %finish : flux calculation
            %(w_bar and V are alive for calculation of FC with/without bump)
                        
            %start : bump analysis
            %extract bump part with spectrum analysis
            
            
            
            
            
            
            
        end % if 18000 skip 
    end % loop for single csv
    clear i temp
    
    
end % loop for all files in data folder
clear t
