% data Ã³¸®3 - daily maximum UHI

% data_hour
% 1: site number
% 2: winddirection
% 3: LCZ
% 4: YYYYMMDDHH
% 5: air temperature
% 6: YYYYMMDD00 (³¯Â¥ - for ÇÇ¹þ)
% 7: UHI (T_AWS - T_±èÆ÷°øÇ×)

% data_kp
% 1: YYYYMMDDHH
% 2: air temperature

% data_UHI
% 1: YYYYMMDD00
% 2: site number
% 3: daily maximum UHI (T_AWS - T_±èÆ÷°øÇ×)
% 4: LCZ at maximum UHI
% 5: HH (time)
% 6: number of all data (24/23 are OK)
% 7: number of good data
% 8: number of gaps

% list_site
% 1: site number
% 2: start number in data_hour
% 3: end number in data_hour
% 4: number of data
% 5: number of days

load('data_hour2.mat');
load('data_kp.mat');

data_hour = sortrows(data_hour,[1 4]);
% 1: site numbber
% 4: date

list_site = unique(data_hour(:,1));
list_site(1,2)=0;
list_site(1,3)=0;
list_site(1,4)=0;
for i = 1:length(list_site)
    temp = find(data_hour(:,1)==list_site(i,1));
    list_site(i,2) = min(temp);
    list_site(i,3) = max(temp);
    list_site(i,4) = list_site(i,3)-list_site(i,2)+1;
end
clear i temp

list_site(1,5)=0;
for i = 1:length(list_site)
    temp = zeros(list_site(i,4),1);
    for j = 1:list_site(i,4)
        temp(j,1) = data_hour(list_site(i,2)-1+j,6);
    end
    clear j
    temp = unique(temp);
    list_site(i,5) = length(temp);
    clear temp
end
clear i 

data_uhi = zeros(sum(list_site(:,5)),8);
n = 0;
for i = 1:length(list_site)
    temp = zeros(list_site(i,4),7); % extract 1-site
    for j = 1:list_site(i,4)
        for k = 1:7
            temp(j,k) = data_hour(list_site(i,2)-1+j,k);
        end
        clear k
    end
    clear j
    temp2 = unique(temp(:,6)); %date
    temp3  =0;
    for j = 1:length(temp2)
        n = n+1;
        data_uhi(n,1) = temp2(j,1); %YYYYMMDD00
        temp3 = find(temp(:,6)==temp2(j,1)); %all positions of the day
        temp3(1,2) = 0; %YYYYMMDDHH
        temp3(1,3) = 0; %air temperature
        temp3(1,4) = 0; %flag for good data
        for k = 1:length(temp3(:,1))
            temp3(k,2) = temp(temp3(k,1),4);
            temp3(k,3) = temp(temp3(k,1),5);
            if (temp3(k,3)<-30)||(temp3(k,3)>45)
                temp3(k,3)= -888.8;
            else
                temp4 = find(data_kp(:,1)==temp3(k,2));
                temp4 = data_kp(temp4,2);
                if (temp4<-30)||(temp4>45)
                    temp3(k,3)= -888.8;
                else
                    temp3(k,4)= 1; %flag
                    temp3(k,3)= temp3(k,3)-temp4; %UHI
                    temp3(k,2)= temp(temp3(k,1),4); %YYYYMMDDHH
                end
            end
        end
        clear k
        if sum(temp3(:,4))>0
            temp4 = find(temp3(:,3)==max(temp3(:,3)));
            data_uhi(n,2) = list_site(i,1);
            data_uhi(n,3) = temp3(temp4(1),3); %max UHI
            data_uhi(n,5) = temp3(temp4(1),2); % HH
            temp4 = find(temp(:,4)==data_uhi(n,5));
            data_uhi(n,4) = temp(temp4,3); % LCZ
            data_uhi(n,6) = length(temp3); % number of all
            data_uhi(n,7) = sum(temp3(:,4)); % number of good
            data_uhi(n,8) = data_uhi(n,6)-data_uhi(n,7); % number of gaps
        end
    end
    clear j
    clear temp temp2 temp3 temp4
end
clear i n


save('data_hour3.mat','data_hour');
save('data_uhi.mat','data_uhi');



