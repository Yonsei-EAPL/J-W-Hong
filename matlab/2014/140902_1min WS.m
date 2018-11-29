%% History
% 140804 Mr. Je-Woo Hong; for BoSeong Tower, CSAT3 only


%% input(*.csv) information

sonic_ang = 90+7.24; %for BS Tower

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

result = zeros(max(num_30min_n(:,1))*30,3);
% 1; n_data (unitless)
% 2; mean wind speed, U (m/s)


%% position
po_u = 5;
po_v = 6;


%% main process
for i = 1:5184000/1200
    result(i,1) = 1200;    
    
    % extract
    temp = zeros(60*20,2);
    for j = 1:1200
        if i == 1
            temp(j,1) = data(j,po_u);
            temp(j,2) = data(j,po_v);
        else
            temp(j,1) = data((i-1)*1200+j,po_u);
            temp(j,2) = data((i-1)*1200+j,po_v);
        end
    end
    clear j
    
    % mean wind-speed, U
    temp_ws = 0;
    temp_n = 0;
    temp_rm = 0;
    for j = 1:1200
        temp_rm = temp_rm+1;
        if abs(temp(temp_rm,1))>20%65.535
            temp(temp_rm,:)=[];
            temp_rm = temp_rm-1;
        elseif abs(temp(temp_rm,2))>20%65.535
            temp(temp_rm,:)=[];
            temp_rm = temp_rm-1;
        else
            temp_ws = temp_ws + (temp(temp_rm,1)^2 + temp(temp_rm,2)^2 )^(0.5);
            temp_n = temp_n + 1;
        end
    end
    temp_ws = temp_ws/temp_n;
    result(i,2) = temp_ws; % result : mean wind-speed, U
    result(i,3) = temp_n;
    clear j temp_ws
end
clear i temp temp_n temp_rm
clear po_u po_v po_w po_Ts 
clear size_n size_var sonic_ang num_30min_n

