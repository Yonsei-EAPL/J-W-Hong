%% History
% 140804 Mr. Je-Woo Hong; for BoSeong Tower, CSAT3 only


%% input(*.csv) information
dataDir = 'E:\EAPL\JW_Observation\보성\csv\140m'; % folder name
% dataDir = 'E:\EAPL\JW_Observation\은평뉴타운\_csv'; % folder name
dataName = dir(fullfile(dataDir, 'CSV_8973.*'))'; 
total_result=cell(length(dataName),2);

wait = waitbar(0,'...Flux Data Process...');


%% position

po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8;


%% input(*.csv) information

sonic_ang = 90+7.24; %for BS Tower

for t=6:length(dataName) %n_st:n_en % 여러개의 CSV 파일을 연속적으로 처리하기 위한 loop

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
            if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==0.05)
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
    % ...
    % 29; cov_wTs (K*m/s)
    % 30; 1200*cov_wTs (W/m2)
    % ...
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



    %% main process
    for i = 1:max(num_30min_n(:,1))
        result(i,1) = num_30min_n(i,2);    
        waitbar((t/length(dataName)),wait,sprintf('%3f',i))        

        % extract
        temp = zeros(num_30min_n(i,2),4);
        for j = 1:num_30min_n(i,2)
            if i == 1
                temp(j,1) = data(j,po_u);
                temp(j,2) = data(j,po_v);
                temp(j,3) = data(j,po_w);
                temp(j,4) = data(j,po_Ts);
            else
                temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
                temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
                temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
                temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
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

        % sensible heat flux ; AmeriFlux
        w_bar = result(i,5); % m/s
        cov_wTs = 0;
        for j = 1:num_30min_n(i,2)
            cov_wTs = cov_wTs + (temp(j,3)-w_bar)*(temp(j,4)-Ts);
        end
        clear j
        cov_wTs = cov_wTs/num_30min_n(i,2);
        result(i,29) = cov_wTs;
        result(i,30) = 1200*cov_wTs;

        clear Ts w_bar cov_wTs
    end
    total_result{t,1}= dataName(t).name;
    total_result{t,2}= result;
    
    clear i
    clear size_n size_var num_30min_n    
end
close(wait);
clear data dataDir dataName num_30min_n result wait t
clear po_u po_v po_w po_Ts sonic_ang

