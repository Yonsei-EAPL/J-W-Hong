%% After Ogive, General Analysis for "1day TOA5 data"
% History
% 2013-06-10, J-W Hong, for Eun-Pyeong and Seoul-Forest, using 1day TOA5
site = 1;    % 1 for EP (CSAT angle = 220), 2 for SF (CSAT angle = 230)

%% Information
Ogive = result;
clear result
if site==1
    sonic_ang = 220;    % for EP
elseif site==2
    sonic_ang = 230;    % for SF
end
clear site
length = 24; % data length in hour
Hz = 10; % sampling rate
avg_period=30;  % averaging period in [min]
var = 12; % number of input variables
po_u = 2; % column number of u
    % po_v = 3;
    % po_w = 4;
    % po_Ts = 5;
    % po_diag_sonic = 6;
    % po_co2 = 7; % in ppm
    % po_h2o = 8; % in g m-3
    % po_diag_irga = 9;
    % po_T = 10;   % irga_T
    % po_P = 11;  % irga_P
    % po_sig_str_co2 = 12;
    % po_sig_str_h2o = 13;
result=zeros(length*(60/avg_period),17); 
    % 1: mean wind-speed
    % 2: mean wind-direction
    % 3: u*
    % 4: std_w/u*
    % 5: mean Ts
    % 6: mean co2 in ppm
    % 7: mean h2o in g m-3
    % 8: max diag_sonic
    % 9: max diag_irga
    %10: min sig_str_co2
    %11: min sig_str_h2o
    %12: mean cell_T
    %13: mean cell_P
    %14: min co2
    %15: max co2
    %16: min h2o
    %17: max h2o

%% Analysis
for i = 1:length*(60/avg_period)

    % extract data
    temp = zeros(avg_period*60*Hz,12);
    for j = 1:avg_period*60*Hz
        for k = 1:var
            temp(j,k) = data((i-1)*avg_period*60*Hz+j,po_u+(k-1));
        end
    end
    clear j k
    
    % mean wind-speed
    temp_ws =0; % for mean wind-speed
    for j = 1:avg_period*60*Hz
        temp_ws = temp_ws + (temp(j,1)^2 + temp(j,2)^2 + temp(j,3)^2)^(0.5);
    end
    temp_ws = temp_ws / (avg_period*60*Hz);
    result(i,1) = temp_ws;
    clear j temp_ws
    
    % mean wind-direction
    %   using u, v component only without rotation
    temp_wd =0; % for mean wind-direction
    for j = 1:avg_period*60*Hz
        if temp(j,1)>0
            if temp(j,2)>0
                temp_wd = temp_wd + 360 - atan(temp(j,2)/temp(j,1))/pi()*180;
            else
                temp_wd = temp_wd + atan((-1*temp(j,2))/temp(j,1))/pi()*180;
            end
        else
            if temp(j,2)>0
                temp_wd = temp_wd + 270 + atan(temp(j,2)/(-1*temp(j,1)))/pi()*180;
            else
                temp_wd = temp_wd + 270 - atan(temp(j,2)/temp(j,1))/pi()*180;
            end
        end
    end
    temp_wd = temp_wd / (avg_period*60*Hz);
    if (temp_wd + sonic_ang) > 360
        temp_wd = temp_wd + sonic_ang - 360;
    else
        temp_wd = temp_wd + sonic_ang;
    end
    result(i,2) = temp_wd;
    clear j temp_wd
    
    % double rotation
    temp_wind = zeros(avg_period*60*Hz,3);  % rotated wind components
    u_bar = mean(temp(:,1)); % for double rotation
    v_bar = mean(temp(:,2));
    w_bar = mean(temp(:,3));
    alpha = atan(v_bar/u_bar);
    beta = atan(w_bar/(sqrt(u_bar^2+v_bar^2)));
    for j = 1:avg_period*60*Hz
        temp_wind(j,1) = cos(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+sin(beta)*temp(j,3);
        temp_wind(j,2) = -sin(alpha)*temp(j,1)+cos(alpha)*temp(j,2);
        temp_wind(j,3) = -sin(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+cos(beta)*temp(j,3);
    end
    clear u_bar v_bar w_bar alpha beta j
    
    % u*
    %   using, u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.5)
    %   ref. Stull, p. 67
    uw = 0; % for mean(u'w')
    vw = 0; % for mean(v'w')
    u_bar = mean(temp_wind(:,1));
    v_bar = mean(temp_wind(:,2));
    w_bar = mean(temp_wind(:,3));
    for j = 1:avg_period*60*Hz
        uw = uw + (temp_wind(j,1)-u_bar)*(temp_wind(j,3)-w_bar);
        vw = vw + (temp_wind(j,2)-v_bar)*(temp_wind(j,3)-w_bar);
    end
    uw = uw/(avg_period*60*Hz);
    vw = vw/(avg_period*60*Hz);
    result(i,3) = (uw^2 + vw^2)^(0.5);
    clear uw vw u_bar v_bar w_bar j
    
    % std_w/u*
    sigma_w = std(temp_wind(:,3));
    result(i,4) = result(i,3)/sigma_w;
    clear sigma_w;
    
    % mean Ts
    result(i,5) = mean(temp(:,4));
    
    % mean co2
    result(i,6) = mean(temp(:,6));
    result(i,14) = min(temp(:,6));
    result(i,15) = max(temp(:,6));
    
    % mean h2o
    result(i,7) = mean(temp(:,7));
    result(i,16) = min(temp(:,7));
    result(i,17) = max(temp(:,7));
    
    % max diag_sonic
    result(i,8) = max(temp(:,5));
    
    % max diag_irga
    result(i,9) = max(temp(:,8));
    
    % min sig_str_co2
    result(i,10) = min(temp(:,11));
    
    % min sig_str_h2o
    result(i,11) = min(temp(:,12));
    
    % mean cell T
    result(i,12) = mean(temp(:,9));
    
    % mean cell P
    result(i,13) = mean(temp(:,10));
    
end
clear temp temp_wind i
clear length var Hz avg_period po_u site sonic_ang