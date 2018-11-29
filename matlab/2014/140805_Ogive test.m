%% History
% 140804 Mr. Je-Woo Hong; for BoSeong Tower, CSAT3 only


%% input(*.csv) information
sonic_ang = 90+7.24; %for BS Tower
SR = 0.05; % sampling rate (sec)

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
        if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==SR)
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
result = zeros(16, max(num_30min_n(:,1))/4);
% 1; cov(wTs) of 1-min averaging
% 2; 2
% 3; 3
% 4; 4
% 5; 5
% 6; 6
% 7; 8
% 8; 10
% 9; 12
% 10; 15
% 11; 20
% 12; 24
% 13; 30
% 14; 40
% 15; 60
% 16; 120


%% position
po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8;


%% running information
period = [1;2;3;4;5;6;8;10;12;15;20;24;30;40;60;120]; % 16 set
n_period = 16;


%% main process
h = waitbar(0,'Please wait...');
for i = 1:max(num_30min_n(:,1))/4
    waitbar(i/(max(num_30min_n(:,1))/4))
    % extract
    temp = zeros(120*60/SR,4);
    for j = 1:120*60/SR
            temp(j,1) = data((i-1)*120*60/SR+j,po_u);
            temp(j,2) = data((i-1)*120*60/SR+j,po_v);
            temp(j,3) = data((i-1)*120*60/SR+j,po_w);
            temp(j,4) = data((i-1)*120*60/SR+j,po_Ts);
    end
    clear j

    % main-run for kinematic sensible heat flux
    for j = 1:n_period
        n_inter = 120/period(j,1);
        for k = 1:n_inter
            % extract
            temp2 = zeros(period(j,1)*60/SR,4);
            for l =1:period(j,1)*60/SR
                temp2(l,1) = temp((k-1)*period(j,1)*60/SR+l,1);
                temp2(l,2) = temp((k-1)*period(j,1)*60/SR+l,2);
                temp2(l,3) = temp((k-1)*period(j,1)*60/SR+l,3);
                temp2(l,4) = temp((k-1)*period(j,1)*60/SR+l,4);
            end
            clear l
            
            % double rotation
            u_bar = mean(temp2(:,1));
            v_bar = mean(temp2(:,2));
            w_bar = mean(temp2(:,3));            
            temp_wind = zeros(period(j,1)*60/SR,3);
            alpha = atan2(v_bar,u_bar);
            beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.5)));
            for l = 1:period(j,1)*60/SR
                temp_wind(l,1) = cos(beta)*(cos(alpha)*temp2(l,1)+sin(alpha)*temp2(l,2))+sin(beta)*temp2(l,3);
                temp_wind(l,2) = -sin(alpha)*temp2(l,1)+cos(alpha)*temp2(l,2);
                temp_wind(l,3) = -sin(beta)*(cos(alpha)*temp2(l,1)+sin(alpha)*temp2(l,2))+cos(beta)*temp2(l,3);
                temp2(l,1) = temp_wind(l,1);
                temp2(l,2) = temp_wind(l,2);
                temp2(l,3) = temp_wind(l,3);
            end
            clear u_bar v_bar alpha beta temp_wind l

            % mean Ts
            Ts = mean(temp2(:,4));

            % kinematic sensible heat flux
            cov_wTs = 0;
            for l = 1:period(j,1)*60/SR
                cov_wTs = cov_wTs + (temp2(l,3)-w_bar)*(temp2(l,4)-Ts);
            end
            clear l
            cov_wTs = cov_wTs/(period(j,1)*60/SR);
            result(j,i) = result(j,i) + cov_wTs;
        end
        result(j,i) = result(j,i)/n_inter;
        clear k
    end
    clear j
end
close(h)
clear i temp h
clear po_u po_v po_w po_Ts 
clear size_n size_var sonic_ang num_30min_n SR
