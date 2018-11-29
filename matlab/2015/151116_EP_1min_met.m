length = 1938;%4913;%2*24*60; %always check this value in each run
sonic_ang = 220+8.06; %for EP NewTown

%% result
result = zeros(length,11);
% 1; year
% 2; DOY
% 3; hhmm
% 4; mean wind speed, U (m/s)
% 5; u_bar (m/s)
% 6; v_bar (m/s)
% 7; w_bar (m/s)
% 8; wind direction (including sonic_angle) (degree)
% 9; standard deviation of wind direction (including sonic_angle) (degree)
% 10; mean Ts (degree C)
% 11; mean Tk (degree C)

%% main process
for i = 1:length
    result(i,1) = data((i-1)*600+300,1);
    result(i,2) = data((i-1)*600+300,2);
    result(i,3) = data((i-1)*600+300,3);

    temp = zeros(600,6);
    for j = 1:600
        temp(j,1) = data((i-1)*600+j,5);
        temp(j,2) = data((i-1)*600+j,6);
        temp(j,3) = data((i-1)*600+j,7);
        temp(j,4) = data((i-1)*600+j,8);
        temp(j,5) = data((i-1)*600+j,11);        
    end
    clear j

    % mean wind-speed, U
    temp_ws = 0;
    for j = 1:600
        temp_ws = temp_ws + (temp(j,1)^2 + temp(j,2)^2 + temp(j,3)^2)^(0.5);
    end
    temp_ws = temp_ws/600;
    result(i,4) = temp_ws; % result : mean wind-speed, U
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
    for j = 1:600
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
    temp_std_wd = (temp_std_wd/600)^(0.5);
    result(i,9) = temp_std_wd; % result : standard deviation of wind-direction
    if mean_wd+sonic_ang>360
        mean_wd = mean_wd + sonic_ang-360;
    else
        mean_wd = mean_wd + sonic_ang;
    end
    result(i,8) = mean_wd; % result : wind direction
    clear j temp_wd temp_std_wd mean_wd 

    result(i,5) = u_bar; % result : mean u
    result(i,6) = v_bar; % result : mean v
    result(i,7) = w_bar; % result : mean w

    % mean Ts
    result(i,10) = mean(temp(:,4)); % result : mean_Ts 
    Ts = result(i,10); % degree C

    % mean h2o
    Xv = mean(temp(:,5)); % result : mean h2o
    Tk = (Ts + 273.15)/(1+0.000321*Xv); % K ; AmeriFlux
    result(i,11) = Tk-273.15; % result : actual temperature in degree C    

end
clear i
clear sonic_ang w_bar v_bar u_bar temp length Xv Ts Tk
