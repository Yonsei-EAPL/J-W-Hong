%% History
% 130714 Mr.Je-Woo Hong, Mr.Keunmin Lee and Prof.Jinkyu Hong

%% using data (EP 2009-2011 comparison)

sonic_ang = 225+8; %for EP NewTown

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

result = zeros(max(num_30min_n(:,1)),6);
% 1; mean wind speed
% 2; wind direction (including sonic_angle)
% 3; u*
% 4; std_w/u*
% 5; mean Ts
% 6; cov_wTs
% 7; mean_u
% 8; std_v
% 9; std_wd
% 10; n_data

%% position

po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8; %for 2010-2012
%po_Ts = 10; %for 2008-2009

%% 

for i = 1:max(num_30min_n(:,1))
    % extract
    temp = zeros(num_30min_n(i,2),4);
    for j = 1:num_30min_n(i,2)
        if i ==1
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
    result(i,2) = mean_wd;
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
    result(i,9) = temp_std_wd;
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
    clear u_bar v_bar w_bar alpha beta j temp_wind

    % u* ; u* = ( mean(u'w')^2 + mean(v'w')^2 )^(0.25) ; Stull, p.67
    uw = 0; % for mean(u'w')
    vw = 0; % for mean(v'w')
    u_bar = mean(temp(:,1));
    result(i,7) = u_bar;
    v_bar = mean(temp(:,2));
    w_bar = mean(temp(:,3));
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

    sigma_v = std(temp(:,2));
    result(i,8) = sigma_v;
    clear sigma_v
    
    % mean Ts
    result(i,5) = mean(temp(:,4));
    
    % sensible heat flux ; AmeriFlux standard
    w_bar = mean(temp(:,3));
    Ts = mean(temp(:,4));
    cov_wTs = 0;
    for j = 1:num_30min_n(i,2)
        cov_wTs = cov_wTs + (temp(j,3)-w_bar)*(temp(j,4)-Ts);
    end
    clear j
    cov_wTs = cov_wTs/num_30min_n(i,2);
    result(i,6) = cov_wTs;
    clear cov_wTs w_bar Ts
    
    result(i,10) = num_30min_n(i,2);
end
clear i temp
clear po_u po_v po_w po_Ts size_n size_var

