%% 2013-08-14 for daily&monthly mean T and accumulated rainfall
% using mod_seoul (after QC/QA test)
% in mod_seoul(n_st, n_data, n_var)
% save_seoul = mod_seoul; % conserv the raw data

[n_st n_data n_var] = size(mod_seoul);
n_days = n_data/24;
n_year = n_days/366; % using 29th Feb.


%% define constants

n_md = [31;29;31;30;31;30;31;31;30;31;30;31];
LCZ = [7;20;15;4;4;1;1;2;10;16;12;11;3;9;8;18;5;12;1;18;3;21;15;19;18;6;4;18;4;14;13;17;4;20];
n_LCZ = 21;
n_amd = zeros(12,1);
for i = 1:12
    for j = 1:i
        n_amd(i,1) = n_amd(i,1)+n_md(j,1);
    end
end
clear i j


%% make bin for result

mean_A_T = zeros(n_st, n_year);
mean_M_T = zeros(n_st,12*n_year);
mean_D_T = zeros(n_st, n_days);

accu_A_R = zeros(n_st, n_year);
accu_M_R = zeros(n_st,12*n_year);
accu_D_R = zeros(n_st, n_days);

UHI_seoul = zeros(n_st, n_days);
mean_M_UHI = zeros(n_st,12*n_year);
mean_M_UHI_LCZ = zeros(n_LCZ,12*n_year);
DTR_seoul = zeros(n_st, n_days);
mean_M_DTR = zeros(n_st,12*n_year);
mean_M_DTR_LCZ = zeros(n_LCZ, 12*n_year);


%% daily process

for i = 1:n_days
    for k = 1:n_st
        temp1=0; % count for T
        for j = 1:24
            if mod_seoul(k,(i-1)*24+j,2)>-90 % bigger than error
                mean_D_T(k,i) = mean_D_T(k,i)+ mod_seoul(k,(i-1)*24+j,2);
                temp1 = temp1+1;
            end
            if mod_seoul(k,(i-1)*24+j,3)>0 % rainfall is positive
                accu_D_R(k,i) = accu_D_R(k,i)+ mod_seoul(k,(i-1)*24+j,3);
            end
        end
        if temp1 == 0
            mean_D_T(k,i) = -999;
        else
            mean_D_T(k,i) = mean_D_T(k,i)/temp1;
        end
    end
end
clear i j k temp1


%% monthly process

temp_m = 0;
temp_nT = zeros(n_st,n_year*12);
temp_nR = zeros(n_st,n_year*12);
for k = 1:n_st
    for j = 1:n_year
        for i = 1:366
            if i==1
                temp_m = 1;
            elseif i > n_amd(temp_m,1)
                temp_m = temp_m + 1;
            end
            if mean_D_T(k,(j-1)*366+i)>-900
                mean_M_T(k,(j-1)*12+temp_m) = mean_M_T(k,(j-1)*12+temp_m)+mean_D_T(k,(j-1)*366+i);
                temp_nT(k,(j-1)*12+temp_m) = temp_nT(k,(j-1)*12+temp_m)+1;
            end
            if accu_D_R(k,(j-1)*366+i)>0
                accu_M_R(k,(j-1)*12+temp_m) = accu_M_R(k,(j-1)*12+temp_m)+accu_D_R(k,(j-1)*366+i);
                temp_nR(k,(j-1)*12+temp_m) = temp_nR(k,(j-1)*12+temp_m)+1;
            end               
        end
    end
    for j = 1:12*n_year
        if temp_nT(k,j)>0
            mean_M_T(k,j) = mean_M_T(k,j)/temp_nT(k,j);
        end
    end
end
n_mdT = temp_nT;
n_mdR = temp_nR;
clear i j k temp_m temp_nT temp_nR


%% annual process
temp_nT = zeros(n_st,n_year);
temp_nR = zeros(n_st,n_year);
for k = 1:n_st
    for i = 1:n_year
        for j = 1:12
            mean_A_T(k,i) = mean_A_T(k,i)+mean_M_T(k,(i-1)*12+j)*n_mdT(k,(i-1)*12+j);
            temp_nT(k,i) = temp_nT(k,i)+n_mdT(k,(i-1)*12+j);
            accu_A_R(k,i) = accu_A_R(k,i)+accu_M_R(k,(i-1)*12+j);
            temp_nR(k,i) = temp_nR(k,i)+n_mdR(k,(i-1)*12+j);
        end
    mean_A_T(k,i) = mean_A_T(k,i)/temp_nT(k,i);
    end
    
end
n_adT = temp_nT;
n_adR = temp_nR;
clear i j k temp_nT temp_nR


%% UHI (urban heat island)
for i = 1:n_st
    for j = 1:n_days
        temp_doy=mod(j,366);
        if temp_doy==0
            temp_doy=366;
        end
        UHI_seoul(i,j) = mod_seoul(i,(j-1)*24+UHI(temp_doy,1),2);
    end
end
clear i j temp_doy
gap_UHI = zeros(n_days,1);
for i = 1:n_days
    for j = 1:34
        if UHI_seoul(j,i)<-900
            gap_UHI(i,1) = gap_UHI(i,1)+1;
        end
    end
end
clear i j

temp_y = 1;
temp_m = 1;
n_mdUHI = zeros(n_st,n_year*12);
for i = 1:n_days
    temp_doy=mod(i,366);
    if temp_doy==0
        temp_doy=366;
    end
    if temp_doy>n_amd(temp_m,1)
        temp_m = temp_m + 1;
    end
    if temp_y<4
        if gap_UHI(i,1)==3
            temp = sum(accu_D_R(:,i));
            if temp==0
                for j = 1:n_st
                    if ((j~=3)&&(j~=23))&&(j~=25)
                        if UHI_seoul(j,i)>-990
                            mean_M_UHI(j,(temp_y-1)*12+temp_m) = mean_M_UHI(j,(temp_y-1)*12+temp_m)+UHI_seoul(j,i);
                            n_mdUHI(j,(temp_y-1)*12+temp_m) = n_mdUHI(j,(temp_y-1)*12+temp_m)+1;
                        end
                    end
                end
            end
            if temp_doy ==366
                temp_y = temp_y + 1;
                temp_m = 1;
            end
        end
    else
        if gap_UHI(i,1)==0
            temp = sum(accu_D_R(:,i));
            if temp==0
                for j = 1:n_st
                    if ((j~=3)&&(j~=23))&&(j~=25)
                        if UHI_seoul(j,i)>-990
                            mean_M_UHI(j,(temp_y-1)*12+temp_m) = mean_M_UHI(j,(temp_y-1)*12+temp_m)+UHI_seoul(j,i);
                            n_mdUHI(j,(temp_y-1)*12+temp_m) = n_mdUHI(j,(temp_y-1)*12+temp_m)+1;
                        end
                    end
                end
            end
            if temp_doy ==366
                temp_y = temp_y + 1;
                temp_m = 1;
            end            
        end
    end
end
for i = 1:n_st
   if ((j~=3)&&(j~=23))&&(j~=25)
        for j = 1:n_year*12
            if n_mdUHI(i,j)>0
                mean_M_UHI(i,j) = mean_M_UHI(i,j)/n_mdUHI(i,j);
            else
                mean_M_UHI(i,j) = -999;
            end
        end
    end
end
clear i j temp temp_m temp_doy temp_y

for i = 1:n_year*12
    temp_n = 0;
    temp_mean =0;
    for j = 1:n_st
        if ((j~=3)&&(j~=23))&&(j~=25)
            if mean_M_UHI(j,i)~=-999
                temp_mean = temp_mean + mean_M_UHI(j,i);
                temp_n = temp_n + 1;
            end
        end
    end
    temp_mean = temp_mean/temp_n;
    temp_nm = zeros(n_LCZ,1);
    for j = 1:n_st
        if ((j~=3)&&(j~=23))&&(j~=25)
            if mean_M_UHI(j,i)~=-999
                mean_M_UHI_LCZ(LCZ(j,1),i) = mean_M_UHI_LCZ(LCZ(j,1),i)+(mean_M_UHI(j,i)-temp_mean);
                temp_nm(LCZ(j,1),1)= temp_nm(LCZ(j,1),1)+1;
            end
        end
    end
    for j = 1:n_LCZ
        if temp_nm(j,1)>0
            mean_M_UHI_LCZ(j,i) = mean_M_UHI_LCZ(j,i)/temp_nm(j,1);
        else
            mean_M_UHI_LCZ(j,i) = -999;
        end
    end
    
end
clear i j temp_mean temp_n temp_nm


%% DTR (diurnal temperature range)

temp_min = -99;
temp_max = -99;
min_T = zeros(n_st,n_days);
max_T = zeros(n_st,n_days);
for i = 1:n_st
    for j = 1:n_days
        min_T(i,j) = -999;
        max_T(i,j) = -999;
    end
end

for i = 1:n_st
    for j = 1:n_days
        for k = 1:24
            if temp_min==-99
                if mod_seoul(i,(j-1)*24+k,2)>-990
                    temp_min = mod_seoul(i,(j-1)*24+k,2);
                    temp_max = mod_seoul(i,(j-1)*24+k,2);
                    min_T(i,j) = mod_seoul(i,(j-1)*24+k,2);
                    max_T(i,j) = mod_seoul(i,(j-1)*24+k,2);
                end
            else
                if (mod_seoul(i,(j-1)*24+k,2)>-990)&&(temp_min>mod_seoul(i,(j-1)*24+k,2))
                    temp_min = mod_seoul(i,(j-1)*24+k,2);
                    min_T(i,j) = mod_seoul(i,(j-1)*24+k,2);
                end
                if (mod_seoul(i,(j-1)*24+k,2)>-990)&&(temp_max<mod_seoul(i,(j-1)*24+k,2))
                    temp_max = mod_seoul(i,(j-1)*24+k,2);
                    max_T(i,j) = mod_seoul(i,(j-1)*24+k,2);
                end                
            end
        end
        if (temp_min~=-99)&&(temp_max~=-99)
            DTR_seoul(i,j) = temp_max - temp_min;
        else
            DTR_seoul(i,j) = -999;
        end
        temp_min = -99;
        temp_max = -99;
    end
end
clear i j k temp_min temp_max 
gap_DTR = zeros(n_days,1);
for i = 1:n_days
    for j = 1:34
        if DTR_seoul(j,i)<-900
            gap_DTR(i,1) = gap_DTR(i,1)+1;
        end
    end
end
clear i j

temp_y = 1;
temp_m = 1;
n_mdDTR = zeros(n_st,n_year*12);
for i = 1:n_days
    temp_doy=mod(i,366);
    if temp_doy==0
        temp_doy=366;
    end
    if temp_doy>n_amd(temp_m,1)
        temp_m = temp_m + 1;
    end
    if temp_y<4
        if gap_DTR(i,1)==3
            temp = sum(accu_D_R(:,i));
            if temp==0
                for j = 1:n_st
                    if ((j~=3)&&(j~=23))&&(j~=25)
                        if DTR_seoul(j,i)~=-999
                            mean_M_DTR(j,(temp_y-1)*12+temp_m) = mean_M_DTR(j,(temp_y-1)*12+temp_m)+DTR_seoul(j,i);
                            n_mdDTR(j,(temp_y-1)*12+temp_m) = n_mdDTR(j,(temp_y-1)*12+temp_m)+1;
                        end
                    end
                end
            end
            if temp_doy ==366
                temp_y = temp_y + 1;
                temp_m = 1;
            end
        end
    else
        if gap_DTR(i,1)==0
            temp = sum(accu_D_R(:,i));
            if temp==0
                for j = 1:n_st
                    if ((j~=3)&&(j~=23))&&(j~=25)
                        if DTR_seoul(j,i)~=-999
                            mean_M_DTR(j,(temp_y-1)*12+temp_m) = mean_M_DTR(j,(temp_y-1)*12+temp_m)+DTR_seoul(j,i);
                            n_mdDTR(j,(temp_y-1)*12+temp_m) = n_mdDTR(j,(temp_y-1)*12+temp_m)+1;
                        end
                    end
                end
            end
            if temp_doy ==366
                temp_y = temp_y + 1;
                temp_m = 1;
            end            
        end
    end
end
for i = 1:n_st
   if ((j~=3)&&(j~=23))&&(j~=25)
        for j = 1:n_year*12
            if n_mdDTR(i,j)>0
                mean_M_DTR(i,j) = mean_M_DTR(i,j)/n_mdDTR(i,j);
            else
                mean_M_DTR(i,j) = -999;
            end
        end
    end
end
clear i j temp temp_m temp_doy temp_y

for i = 1:n_year*12
    temp_n = 0;
    temp_mean =0;
    for j = 1:n_st
        if ((j~=3)&&(j~=23))&&(j~=25)
            if mean_M_DTR(j,i)~=-999
                temp_mean = temp_mean + mean_M_DTR(j,i);
                temp_n = temp_n + 1;
            end
        end
    end
    temp_mean = temp_mean/temp_n;
    temp_nm = zeros(n_LCZ,1);
    for j = 1:n_st
        if ((j~=3)&&(j~=23))&&(j~=25)
            if mean_M_DTR(j,i)~=-999
                mean_M_DTR_LCZ(LCZ(j,1),i) = mean_M_DTR_LCZ(LCZ(j,1),i)+(mean_M_DTR(j,i)-temp_mean);
                temp_nm(LCZ(j,1),1)= temp_nm(LCZ(j,1),1)+1;
            end
        end
    end
    for j = 1:n_LCZ
        if temp_nm(j,1)>0
            mean_M_DTR_LCZ(j,i) = mean_M_DTR_LCZ(j,i)/temp_nm(j,1);
        else
            mean_M_DTR_LCZ(j,i) = -999;
        end
    end
    
end
clear i j temp_mean temp_n temp_nm

gap_min_T = zeros(n_days,1);
for i = 1:n_days
    for j = 1:n_st
        if min_T(j,i)<-900
            gap_min_T(i,1) = gap_min_T(i,1)+1;
        end
    end
end
clear i j
gap_max_T = zeros(n_days,1);
for i = 1:n_days
    for j = 1:n_st
        if max_T(j,i)<-900
            gap_max_T(i,1) = gap_max_T(i,1)+1;
        end
    end
end
clear i j

mean_min_T = zeros(n_st,1);
mean_max_T = zeros(n_st,1);
temp = 0;
for i = 1:n_days
    if i <= 366*3
        if gap_min_T(i,1)==3
            temp = temp+1;
            for j = 1:n_st
                mean_min_T(j,1) = mean_min_T(j,1)+ min_T(j,i);
            end
        end
    else
        if gap_min_T(i,1)==0
            temp = temp+1;
            for j = 1:n_st
                mean_min_T(j,1) = mean_min_T(j,1)+ min_T(j,i);
            end
        end        
    end
end
for i = 1:n_st
    mean_min_T(i,1) = mean_min_T(i,1)/temp;
end
clear i j temp
temp = 0;
for i = 1:n_days
    if i <= 366*3
        if gap_max_T(i,1)==3
            temp = temp+1;
            for j = 1:n_st
                mean_max_T(j,1) = mean_max_T(j,1)+ max_T(j,i);
            end
        end
    else
        if gap_max_T(i,1)==0
            temp = temp+1;
            for j = 1:n_st
                mean_max_T(j,1) = mean_max_T(j,1)+ max_T(j,i);
            end
        end        
    end
end
for i = 1:n_st
    mean_max_T(i,1) = mean_max_T(i,1)/temp;
end
clear i j temp


% for i = 1:n_st
%     temp1 = 0;
%     temp2 = 0;
%     temp3 = 0;
%     for j = 1:n_year
%         for k = 1:366
%             if k ==1
%                 temp1 = 1;
%             elseif k ==366
%                 temp = sum(accu_D_R(:,(j-1)*366+k));
%                 if temp == 0
%                     if (DTR_seoul(i,(j-1)*366+k)>0)&&(DTR_seoul(i,(j-1)*366+k)<30)
%                         temp2 = temp2 + 1;
%                         temp3 = temp3 + DTR_seoul(i,(j-1)*366+k);
%                     end
%                 end
%                 mean_M_DTR(i,(j-1)*12+temp1) = temp3/temp2;
%                 if temp2 == 0
%                     mean_M_DTR(i,(j-1)*12+temp1) = -999;
%                 end
%                 temp2 = 0;
%                 temp3 = 0;
%             else
%                 if n_amd(temp1,1)<k
%                     mean_M_DTR(i,(j-1)*12+temp1) = temp3/temp2;
%                     if temp2 == 0
%                         mean_M_DTR(i,(j-1)*12+temp1) = -999;
%                     end
%                     temp1 = temp1 +1;
%                     temp2 = 0;
%                     temp3 = 0;
%                 end
%             end
%             temp = sum(accu_D_R(:,(j-1)*366+k));
%             if temp == 0
%                 if (DTR_seoul(i,(j-1)*366+k)>0)&&(DTR_seoul(i,(j-1)*366+k)<30)
%                     temp2 = temp2 + 1;
%                     temp3 = temp3 + DTR_seoul(i,(j-1)*366+k);
%                 end
%             end
%         end
%     end
% end
% clear i j k temp1 temp2 temp3 temp
% 
% temp1 = 0;
% temp2 = zeros(n_LCZ,n_year*12);
% for i = 1:n_year
%     for j = 1:366
%         if j ==1
%             temp1 = 1;
%         else
%             if n_amd(temp1,1)<j
%                 temp1 = temp1 +1;
%             end
%         end
%         temp = sum(accu_D_R(:,(i-1)*366+j));
%         if temp == 0    
%             for k = 1:n_st
%                 if (DTR_seoul(k,(i-1)*366+j)>0)&&(DTR_seoul(k,(i-1)*366+j)<30)
%                     mean_M_DTR_LCZ(LCZ(k,1),(i-1)*12+temp1) = mean_M_DTR_LCZ(LCZ(k,1),(i-1)*12+temp1)+ DTR_seoul(k,(i-1)*366+j);
%                     temp2(LCZ(k,1),(i-1)*12+temp1) = temp2(LCZ(k,1),(i-1)*12+temp1) +1;
%                 end
%             end
%         end
%     end
% end
% for i = 1:n_year*12
%     for j = 1:n_LCZ
%         mean_M_DTR_LCZ(j,i) = mean_M_DTR_LCZ(j,i)/temp2(j,i);
%         if temp2(j,i)==0
%             mean_M_DTR_LCZ(j,i) = -999;
%         end
%     end
% end
% clear i j k temp1 temp2


%% plot
% 
% for i = 1:24
%     figure(1)
%     plot(accu_D_R(i,:),'x')
%     hold on
%     figure(2)
%     plot(mean_D_T(i,:),'x')
%     hold on
%     figure(3)
%     plot(mean_M_T(i,:),'x')
%     hold on
%     figure(4)
%     plot(accu_M_R(i,:),'x')
%     hold on
% end
% clear i

% figure(1)
% plot(mean_M_UHI(13,:),'x')
% hold on
% plot(mean_M_UHI(21,:),'x')
% plot(mean_M_UHI(16,:),'xr')
% plot(mean_M_UHI(20,:),'xr')
% 
% figure(2)
% plot(DTR_seoul(13,:),'x')
% hold on
% plot(DTR_seoul(21,:),'x')
% plot(DTR_seoul(16,:),'xr')
% plot(DTR_seoul(20,:),'xr')

