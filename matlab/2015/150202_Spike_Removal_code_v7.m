tic
%% Description
% L0 -> raw data from data logger
% L1 -> spike removed data
% L2 -> coordinate rotated and processed half-hourly data 



%% Define Constants and Physical Ranges for Spike Detection

thrsh_sonic = 4.0;       % threshold for removing spike using median
thrsh_irga = 2.0;       % threshold for removing spike using median
d_range = 10000;     % processing range as the # of time (10 -> 1s) 

%  SONIC ------------------------------------------------------------------
f_u = 15; % limit of sonic (m/s)
f_v = 15; % limit of sonic (m/s)
f_w = 8.192; % limit of sonic (m/s)
f_Ts_upp = 45.0;   % upper limit of sonic air temperature (C)
f_Ts_low = -15.0;   % upper limit of sonic air temperature (C)

%  IRGA -------------------------------------------------------------------

f_CO2_upp = 700.0; % upper limit of CO2 concent (umol/mol) 
f_CO2_low = 350.0; % lower limit of CO2 concent (umol/mol) 

f_H2O_upp = 40.0;  % upper limit of H2O concent (mmol/mol) 
f_H2O_low = 0;  % lower limit of H2O concent (mmol/mol) 




%% Define Constants for Data Process

sonic_ang = 230+8.04; %for SF sonic_ang = 220+8.06; %for EP NewTown

% Data position
po_u = 5; 
po_v = 6; 
po_w = 7;
po_Ts = 8;
po_diag_son = 9;
po_CO2 = 10;
po_H2O = 11;
po_diag_irga = 12;
%po_cell_tmpr = 13;
po_cell_prs = 14;
%po_CO2_str = 15;
%po_H2O_str = 16;
po_diff_prs = 17; 
 


%% L0 : Load CSV Data (Raw Data)

dataDir = 'E:\EAPL_keun\_EAPL_Observation\서울숲\CSV\new'; % folder name
dataName = dir(fullfile(dataDir, 'CSV_5518.ts_data_*'))'; %for SF
%dataName = dir(fullfile(dataDir, 'CSV_3302.ts_data_*'))'; %for EP NewTown

total_result=cell(length(dataName),2);


for t=1:length(dataName) % 여러개의 CSV 파일을 연속적으로 처리하기 위한 loop

    L0_data = importdata(fullfile(dataDir, dataName(t).name));  
    t2 = length(L0_data);
    
    
    
%% L0 -> L1 : Spike Removal    

    % Load L0 data into 'temp'
    temp = zeros(t2,8);
    temp(:,1) = L0_data(:,po_u); % u
    temp(:,2) = L0_data(:,po_v); % v
    temp(:,3) = L0_data(:,po_w); % w
    temp(:,4) = L0_data(:,po_Ts); % Ts
    temp(:,5) = L0_data(:,po_CO2); % CO2
    temp(:,6) = L0_data(:,po_H2O); % H2O
    temp(:,7) = L0_data(:,po_cell_prs); % cell_prs
    temp(:,8) = L0_data(:,po_diff_prs); % diff_prs
    
    diag = zeros(t2,2);
    diag(:,1) = L0_data(:,po_diag_son); % diag_son
    diag(:,2) = L0_data(:,po_diag_irga); % diag_irga
        
    %clear L0_data
    
                    
    % Calculation diferences for spike removal
    d2 = zeros(t2,6);
        for k = 1:6
            for j = 2:t2-1
                d2(j,k) = (temp(j,k)-temp(j-1,k)) - (temp(j+1,k)-temp(j,k));
            end
            d2(1,k) = d2(2,k); d2(t2,k)=d2(t2-1,k); % 양 끝값 처리
        end
        
        clear j k    
    
	% 각 30분 데이터를 n_range 설정에 따라 b개로 나누어 Spike Detection하기 위한 작업
	num_bin = zeros(ceil(t2/d_range),3); 
	for x = 1:t2
        if (x==1)
            b = 1;
            num_bin(b,1) = b;
            num_bin(b,2) = 1;
            num_bin(b,3) = 1;
        else
            if (num_bin(b,2)==d_range)
                b = b+1;
                num_bin(b,1) = b;
                num_bin(b,2) = 1;
                num_bin(b,3) = num_bin(b-1,3)+1;
            else
                num_bin(b,2) = num_bin(b,2)+1;
                num_bin(b,3) = num_bin(b,3)+1;
            end
        end 
    end
        
	clear x

	%------------------------------------------------------------
	%                    Main Processing
	%------------------------------------------------------------
	% Physical Range Detection Processing
	flag = zeros(t2,6);    % Spike flag
	NaN_num = zeros(1,6);  % Spike number
    
	for x = 1:b
        b2 = num_bin(x,2); % 각 loop당 데이터 개수
        b3 = num_bin(x,3); % 누적 데이터 개수
        b4 = b3-b2+1;      % 각 loop당 데이터 시작 시점

        for j = b4:b3
                % u
                if(abs(temp(j,1))>=f_u);   flag(j,1) = 1; end
                % v
                if(abs(temp(j,2))>=f_v);   flag(j,2) = 1; end
                % w
                if(abs(temp(j,3))>=f_w);   flag(j,3) = 1; end
                % Ts
                if(temp(j,4)>=f_Ts_upp);   flag(j,4) = 1; end
                if(temp(j,4)<=f_Ts_low);   flag(j,4) = 1; end
                % CO2
                if(temp(j,5)>=f_CO2_upp);  flag(j,5) = 1; end 
                if(temp(j,5)<=f_CO2_low);  flag(j,5) = 1; end                     
                % H2O
                if(temp(j,6)>=f_H2O_upp);  flag(j,6) = 1; end                         
                if(temp(j,6)<=f_H2O_low);  flag(j,6) = 1; end
                % diff_prs
                if(temp(j,8)<-7);        flag(j,5:6) = 1; end                      
                % diag_son
                if(diag(j,1)~=0);        flag(j,1:4) = 1; end
                % diag_irga
                if(diag(j,2)~=0);        flag(j,5:6) = 1; end                      
        end
        
        clear j
        
        dd = zeros(b2,6); ddd = zeros(b2,6);
            for k = 1:6
                dd(:,k) = d2(b4:b3,k);
            end
                
        Md = nanmedian(dd);
            for y = 1:b2
                ddd(y,1:6) = dd(y,1:6)-Md(1,:);
            end
        MAD = nanmedian(abs(ddd));
        cr_sonic = thrsh_sonic/0.6745*MAD;
        cr_irga = thrsh_irga/0.6745*MAD;        

        clear k y
        
        % Spike Detection Processing
        %for k = 1:4  % sonic
        %    for j = b4:b3
        %        if (d2(j,k)<Md(1,k)-cr_sonic(1,k))||(d2(j,k)>Md(1,k)+cr_sonic(1,k))
        %            flag(j,k) = 1;
        %        end       
        %    end
        %end
    
        for k = 5  % irga_CO2
            for j = b4:b3
                if (d2(j,k)<Md(1,k)-cr_irga(1,k))||(d2(j,k)>Md(1,k)+cr_irga(1,k))
                    flag(j,k) = 1;
                end       
            end
        end
        
    	clear dd ddd Md MAD cr_sonic cr_irga j k
        clear b2 b3 b4 
        
    end
    
    
    clear x diag
        
    % Spike Removal Processing
	for  k = 1:4
        for j = 1:t2
            if (flag(j,k)==1)
                temp(j,1:8) = NaN; % Sonic의 문제일땐 전부 NaN
                NaN_num(1,k) = NaN_num(1,k)+1;
            end
        end
    end
    
    clear j k 
    
	for  k = 5:6
        for j = 1:t2
            if (flag(j,k)==1)
                temp(j,5:8) = NaN; % IRGA의 문제일땐 IRGA만 NaN
                NaN_num(1,k) = NaN_num(1,k)+1;
            end
        end
    end
    
    clear j k 
       
    %L1_data = cell(1,2);   
    %L1_data{1,1} = temp;
    %L1_data{1,2} = num_bin;
    %L1_data{1,3} = flag;
    
    L1_data = temp;
    
    clear b b2 b3 b4 num_bin
	clear x d2 temp flag NaN_num


    
    
%% L1 -> L2 : Coordinate Rotation, Correction, Making Half-hourly Data

    % 각 CSV 파일을 n개의 30분 단위의 데이터로 나누어 처리하기 위한 작업
    num_30min = zeros(1,3);
    for i = 1:t2;
        if i==1
            n = 1;
            num_30min(n,1) = n;
            num_30min(n,2) = 1;
            num_30min(n,3) = 1;
        else
            if ((mod(L0_data(i,3),100)==0)||(mod(L0_data(i,3),100)==30))&&(mod(L0_data(i,4),1)==0.1)
                n = n+1;
                num_30min(n,1) = n;
                num_30min(n,2) = 1;
                num_30min(n,3) = num_30min(n-1,3)+1;
            else
                num_30min(n,2) = num_30min(n,2)+1;
                num_30min(n,3) = num_30min(n,3)+1;
            end
        end 
    end
    
    clear i L0_data
    
    
    
    result = zeros(n,65);
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
    % 53; n_data (unitless)
    % 54; mean wind speed, U (m/s)
    % 55; u_bar (m/s)
    % 56; v_bar (m/s)
    % 57; w_bar (m/s)
    % 58; wind direction (including sonic_angle) (degree)
    % 59; std_wd (degree)
    % 60; u* (m/s)
    % 61; cov_uw (m2/s2)
    % 62; alpha (1st, rotation)
    % 63; beta (2nd, rotation)
    % 64; mean Ts (degree C)
    % 65; cov_wTs (K*m/s)
    
    
    % constant
    R = 8.3143*10^(-6); % kPa*m3/K/mmol
    Mc = 12; % mgC/mmol
    Md = 0.029; % g/mmol
    Mv = 0.018015; % g/mmol
    
    
    
    %------------------------------------------------------------
    %                    Main Processing
    %------------------------------------------------------------     
    for i = 1:n % n개의 30분 데이터를 연속적으로 처리하기 위한 loop
        n2 = num_30min(i,2); % 각 loop 데이터 개수
        n3 = num_30min(i,3); % 누적 데이터 개수
        n4 = n3-n2+1;        % 각 loop당 데이터 시작 시점
     
        % Load L1 data into 'temp'
        temp = zeros(n2,8);
        temp(:,1) = L1_data(n4:n3,1); % u
        temp(:,2) = L1_data(n4:n3,2); % v
        temp(:,3) = L1_data(n4:n3,3); % w
        temp(:,4) = L1_data(n4:n3,4); % Ts
        temp(:,5) = L1_data(n4:n3,5); % CO2
        temp(:,6) = L1_data(n4:n3,6); % H2O
        temp(:,7) = L1_data(n4:n3,7); % cell_prs
        temp(:,8) = L1_data(n4:n3,8); % diff_prs
            
                
        % Nan이 아닌 갯수 체크
        temp_NaN = sum(~isnan(temp));
        
        % Sonic값만을 이용한 플럭스처리(cov_wTs)
        if temp_NaN(1,1)~=0 
            n_data_son = temp_NaN(1,1);  
            [x,y] = sort(temp(:,1),'ascend'); temp = temp(y,:); % u 순으로 정렬
            temp_son = temp(1:n_data_son,1:4);
            result(i,53) = n_data_son;                        % result : n_data_son
        
            clear x y
            
            % mean wind-speed, U
            temp_son_ws = 0;
            for j = 1:n_data_son
                temp_son_ws = temp_son_ws + (temp_son(j,1)^2 + temp_son(j,2)^2 + temp_son(j,3)^2)^(0.5);
            end
            temp_son_ws = temp_son_ws/n_data_son;
        
            result(i,54) = temp_son_ws;                       % result : mean wind-speed, U
    
            clear j temp_son_ws
    
        
            % mean wind-direction
            u_bar = mean(temp_son(:,1));
            v_bar = mean(temp_son(:,2));
            w_bar = mean(temp_son(:,3));
            if u_bar>0
                if v_bar>0
                    temp_son_wd = 360 - atan(v_bar/u_bar)/pi()*180;
                else
                    temp_son_wd = atan((-1*v_bar)/u_bar)/pi()*180;
                end
            else
                if v_bar>0
                    temp_son_wd = 180 + atan(v_bar/(-1*u_bar))/pi()*180;
                else
                    temp_son_wd = 180 - atan(v_bar/u_bar)/pi()*180;
                end
            end
            mean_wd = temp_son_wd;

            temp_son_std_wd = 0;
            for j = 1:n_data_son
                if temp_son(j,1)>0
                    if temp_son(j,2)>0
                        temp_son_wd = 360 - atan(temp_son(j,2)/temp_son(j,1))/pi()*180;
                    else
                        temp_son_wd = atan((-1*temp_son(j,2))/temp_son(j,1))/pi()*180;
                    end
                else
                    if temp_son(j,2)>0
                        temp_son_wd = 180 + atan(temp_son(j,2)/(-1*temp_son(j,1)))/pi()*180;
                    else
                        temp_son_wd = 180 - atan(temp_son(j,2)/temp_son(j,1))/pi()*180;
                    end
                end
    
                if abs(mean_wd-temp_son_wd)>180
                    if (mean_wd-temp_son_wd)<0
                        temp_son_std_wd = temp_son_std_wd + (mean_wd-temp_son_wd+360)^2;
                    else
                        temp_son_std_wd = temp_son_std_wd + (360 - mean_wd-temp_son_wd)^2;
                    end
                else
                    temp_son_std_wd = temp_son_std_wd + (mean_wd-temp_son_wd)^2;
                end
            end
            temp_son_std_wd = (temp_son_std_wd/n_data_son)^(0.5);
        
            if mean_wd+sonic_ang>360
                mean_wd = mean_wd + sonic_ang-360;
            else
                mean_wd = mean_wd + sonic_ang;
            end
            result(i,58) = mean_wd;                       % result : wind direction
            result(i,59) = temp_son_std_wd;                   % result : standard deviation of wind-direction
        
            clear j temp_son_wd temp_son_std_wd mean_wd 

    
            % double rotation
            temp_son_wind = zeros(n_data_son,3);
            alpha = atan2(v_bar,u_bar);
            beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.5)));
            for j = 1:n_data_son
                temp_son_wind(j,1) = cos(beta)*(cos(alpha)*temp_son(j,1)+sin(alpha)*temp_son(j,2))+sin(beta)*temp_son(j,3);
                temp_son_wind(j,2) = -sin(alpha)*temp_son(j,1)+cos(alpha)*temp_son(j,2);
                temp_son_wind(j,3) = -sin(beta)*(cos(alpha)*temp_son(j,1)+sin(alpha)*temp_son(j,2))+cos(beta)*temp_son(j,3);
                temp_son(j,1) = temp_son_wind(j,1);
                temp_son(j,2) = temp_son_wind(j,2);
                temp_son(j,3) = temp_son_wind(j,3);
            end
            result(i,62) = alpha;                        % result : 1st rotation angle
            result(i,63) = beta;                         % result : 2nd rotation angle
            
            clear j u_bar v_bar w_bar alpha beta temp_son_wind

        
            % friction velocity
            % u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.25)    % Stull, p.67
            u_bar = mean(temp_son(:,1)); 
            v_bar = mean(temp_son(:,2)); 
            w_bar = mean(temp_son(:,3)); 
            uw = 0; % for mean(u'w')
            vw = 0; % for mean(v'w')
            for j = 1:n_data_son
                uw = uw + (temp_son(j,1)-u_bar)*(temp_son(j,3)-w_bar);
                vw = vw + (temp_son(j,2)-v_bar)*(temp_son(j,3)-w_bar);
            end
            uw = uw/n_data_son;
            vw = vw/n_data_son;
            result(i,55) = u_bar;                         % result : mean u
            result(i,56) = v_bar;                         % result : mean v
            result(i,57) = w_bar;                         % result : mean w
            result(i,60) = (uw^2 + vw^2)^(0.25);          % result : u_star
            result(i,61) = uw;                            % result : cov_uw
                   
            clear j uw vw u_bar v_bar w_bar            
                       
            % mean Ts
            Ts = mean(temp_son(:,4)); % degree C
            result(i,64) = Ts;                           % result : mean_Ts 
            
            % Buoyancy Flux
            w_bar = result(i,5); % m/s
            cov_wTs = 0;
            for j = 1:n_data_son
                cov_wTs = cov_wTs + (temp_son(j,3)-w_bar)*(temp_son(j,4)-Ts);
            end
            clear j
            cov_wTs = cov_wTs/n_data_son;
                        
            result(i,65) = cov_wTs;                      % result : cov_wTs
            
            clear Ts w_bar cov_wTs n_data_son temp_son
        
            
            
            
            
            
            % CO2 필터 후 플럭스처리
            if temp_NaN(1,5)~=0            
            n_data = temp_NaN(1,5);  
            [x,y] = sort(temp(:,5),'ascend'); temp = temp(y,:); % CO2 농도순으로 정렬
            temp = temp(1:n_data,:);
            result(i,1) = n_data;                        % result : n_data
        
            clear x y
            
       
            % mean wind-speed, U
            temp_ws = 0;
            for j = 1:n_data
                temp_ws = temp_ws + (temp(j,1)^2 + temp(j,2)^2 + temp(j,3)^2)^(0.5);
            end
            temp_ws = temp_ws/n_data;
        
            result(i,2) = temp_ws;                       % result : mean wind-speed, U
    
            clear j temp_ws
    
        
            % mean wind-direction
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

            temp_std_wd = 0;
            for j = 1:n_data
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
            temp_std_wd = (temp_std_wd/n_data)^(0.5);
        
            if mean_wd+sonic_ang>360
                mean_wd = mean_wd + sonic_ang-360;
            else
                mean_wd = mean_wd + sonic_ang;
            end
            result(i,6) = mean_wd;                       % result : wind direction
            result(i,7) = temp_std_wd;                   % result : standard deviation of wind-direction
        
            clear j temp_wd temp_std_wd mean_wd 

    
            % double rotation
            temp_wind = zeros(n_data,3);
            alpha = atan2(v_bar,u_bar);
            beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.5)));
            for j = 1:n_data
                temp_wind(j,1) = cos(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+sin(beta)*temp(j,3);
                temp_wind(j,2) = -sin(alpha)*temp(j,1)+cos(alpha)*temp(j,2);
                temp_wind(j,3) = -sin(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+cos(beta)*temp(j,3);
                temp(j,1) = temp_wind(j,1);
                temp(j,2) = temp_wind(j,2);
                temp(j,3) = temp_wind(j,3);
            end
            result(i,10) = alpha;                        % result : 1st rotation angle
            result(i,11) = beta;                         % result : 2nd rotation angle
            result(i,35) = std(temp(:,1));               % result : standard deviation of u
            result(i,36) = std(temp(:,2));               % result : standard deviation of v
            result(i,37) = std(temp(:,3));               % result : standard deviation of w
            result(i,45) = skewness(temp(:,1));          % result : skewness of u
            result(i,46) = skewness(temp(:,2));          % result : skewness of v
            result(i,47) = skewness(temp(:,3));          % result : skewness of w
            result(i,49) = kurtosis(temp(:,1));          % result : kurtosis of u
            result(i,50) = kurtosis(temp(:,2));          % result : kurtosis of v
            result(i,51) = kurtosis(temp(:,3));          % result : kurtosis of w
    
            clear j u_bar v_bar w_bar alpha beta temp_wind

        
            % friction velocity
            % u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.25)    % Stull, p.67
            u_bar = mean(temp(:,1)); 
            v_bar = mean(temp(:,2)); 
            w_bar = mean(temp(:,3)); 
            uw = 0; % for mean(u'w')
            vw = 0; % for mean(v'w')
            for j = 1:n_data
                uw = uw + (temp(j,1)-u_bar)*(temp(j,3)-w_bar);
                vw = vw + (temp(j,2)-v_bar)*(temp(j,3)-w_bar);
            end
            uw = uw/n_data;
            vw = vw/n_data;
            result(i,3) = u_bar;                         % result : mean u
            result(i,4) = v_bar;                         % result : mean v
            result(i,5) = w_bar;                         % result : mean w
            result(i,8) = (uw^2 + vw^2)^(0.25);          % result : u_star
            result(i,9) = uw;                            % result : cov_uw
            result(i,44) = uw/result(i,35)/result(i,37); % result : r_uw
        
            clear j uw vw u_bar v_bar w_bar
    
    
            % std_u,v,w/u*,U
            result(i,38) = result(i,35)/result(i,8);     % result : std_u/u*
            result(i,39) = result(i,36)/result(i,8);     % result : std_v/u*
            result(i,40) = result(i,37)/result(i,8);     % result : std_w/u*
            result(i,41) = result(i,35)/result(i,2);     % result : std_u/U
            result(i,42) = result(i,36)/result(i,2);     % result : std_v/U
            result(i,43) = result(i,37)/result(i,2);     % result : std_w/U
    
        
            % mean Ts
            Ts = mean(temp(:,4)); % degree C
            result(i,12) = Ts;                           % result : mean_Ts 
            result(i,48) = skewness(temp(:,4));          % result : skewness of Ts
            result(i,52) = kurtosis(temp(:,4));          % result : kurtosis of Ts
  
        
            % mean co2
            Xc = mean(temp(:,5)); % umol/mol (dry air mole fraction)
            result(i,16) = Xc;                           % result : mean co2 
            result(i,17) = min(temp(:,5));               % result : min co2
            result(i,18) = max(temp(:,5));               % result : max co2
    
    
            % mean h2o
            Xv = mean(temp(:,6)); % mmol/mol (dry air mole fraction)
            Tk = (Ts + 273.15)/(1+0.000321*Xv); % K ; AmeriFlux
            result(i,19) = Xv;                           % result : mean h2o
            result(i,20) = min(temp(:,6));               % result : min h2o
            result(i,21) = max(temp(:,6));               % result : max h2o
            result(i,13) = Tk-273.15;                    % result : actual temperature in degree C    
    
    
            % mean cell P
            result(i,14) = mean(temp(:,7));              % result : cell pressure in kPa
    
        
            % mean ambient P
            amb_prs = (mean(temp(:,7))-mean(temp(:,8))); % kPa    
            result(i,15)= amb_prs;                       % result : ambient pressure in kPa
    
    
            % vapor pressure
            e = Xv*amb_prs/(1000+Xv); % kPa
            es = 6.1365*exp(17.502*(Tk-273.15)/(240.97+(Tk-273.15)))/10; % kPa    
            RH = e/es*100;         
            result(i,22) = e;                            % result : vapor pressure in kPa
            result(i,23) = es;                           % result : saturation vapor pressure in kPa
            result(i,24) = RH;                           % result : relative humidity in percent
    
        
            % density
            rho_v = Xv*amb_prs*Mv/(R*Tk*(1000+Xv)); % g/m3
            rho_c = Xc*Mc*10^(-6)*(amb_prs/(R*Tk) - rho_v/Mv); % mgC/m3
            rho_d = (amb_prs-e)*Md/(R*Tk); % g/m3
            result(i,25) = rho_d;                        % result : rho_d
            result(i,26) = rho_v;                        % result : rho_v
            result(i,27) = rho_c;                        % result : rho_c
    
    
            % specific heat capacity
            Cpd = 1005+((Tk-250.03)^2/3364); % J/K/kg ; AmeriFlux
            Cpm = 1859 + 0.13*RH + (Tk-273.15)*(0.193+0.00569*RH) + (Tk-273.15)^2*(0.001+0.000005*RH); % J/K/kg ; AmeriFlux
            rhoCp = Cpm*(rho_v/1000) + Cpd*(rho_d/1000); % specific heat capacity (J/K/kg)
            result(i,28) = rhoCp; % result : rho_Cp (J/K/kg)
    
        
            % sensible heat flux ; AmeriFlux
            w_bar = result(i,5); % m/s
            cov_wTs = 0;
            cov_wq = 0;
            for j = 1:n_data
                cov_wTs = cov_wTs + (temp(j,3)-w_bar)*(temp(j,4)-Ts);
                cov_wq = cov_wq + (temp(j,3)-w_bar)*(temp(j,6)-Xv); % mmol/mol *m/s
            end
            clear j
            cov_wTs = cov_wTs/n_data;
            cov_wq = cov_wq/n_data;
            Qh = rhoCp * (cov_wTs - (0.000321*Tk*cov_wq)); % AmeriFlux
            result(i,29) = cov_wTs;                      % result : cov_wTs
            result(i,30) = rhoCp*cov_wTs;                % result : rho_Cp*cov_wTs
            result(i,31) = Qh;                           % result : Qh
    
    
            % latent heat flux ; AmeriFlux
            lambda = 2500.8 - 2.3668*Ts; % J/g
            V = 0.082*Tk*10^(-3)/(amb_prs/101.325); % m3/mol
            Qe = lambda*Mv*cov_wq/V; 
            result(i,32) = lambda;                       % result : lambda
            result(i,33) = Qe;                           % result : Qe
    
        
            % co2 flux ; AmeriFlux
            co2 = result(i,16); % mean co2 (umol/mol)
            cov_wc = 0;
            for j = 1:n_data
                cov_wc = cov_wc + (temp(j,3)-w_bar)*(temp(j,5)-co2);
            end
            clear j 
            cov_wc = cov_wc/n_data; % umol/mol * m/s
            Fc = cov_wc/V; % umol/m2/s 
            result(i,34) = Fc;                           % result : Fc   

            clear Xv Xc amb_prs
            clear Ts q Tk cell_prs e Pd rho_v rho_d rho_c Cpd es RH Cpm rhoCp w_bar cov_wTs cov_wq Qh
            clear lambda V Qe 
            clear co2 cov_wc Fc
            end
            clear n2 n3 n4 n_data temp temp_NaN
        end                    
    end
    
    clear i n
    clear R Mc Md Mv
            
   %%
    %L2_data = cell(1,2);   
    %L2_data{1,1} = result;
    %L2_data{1,2} = num_30min;
    %clear result num_30min
        
    csvwrite (dataName(t).name, L1_data)
    clear L1_data 
    
    total_result{t,1}= dataName(t).name;
    total_result{t,2}= result;
    %total_result{t,3}= num_30min;
    

    clear result num_30min t2

    %clear L2_data

end

clear sonic_ang t dataName dataDir
clear po_u po_v po_w po_Ts po_CO2 po_H2O po_cell_prs po_diff_prs ...
        po_diag_son po_diag_irga
clear d_range thrsh_irga thrsh_sonic
clear f_CO2_low f_CO2_upp f_H2O_low f_H2O_upp f_Ts_low f_Ts_upp ...
        f_u f_v f_w
    
toc   