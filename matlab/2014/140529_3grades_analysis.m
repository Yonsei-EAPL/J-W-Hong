%%
% time from excel (data2.mat)
i_stn = [2;1;1;1;3;3;3;2;1;1;1;1;3;2;3;1;1;1;1;2;1;2;1;2;3;2;3;1;1;1;1;1;3;3;2;];
n_data = 3988320;
n_data1 = 113952;
n_lcz = 3;
n_stn = [18;8;9];

for i = 1:n_data1
    if (time(i,1) <3)||(time(i,1) ==12)
        time(i,4) = 4;
    elseif (time(i,1)>2)&&(time(i,1) <6)
        time(i,4) = 1;
    elseif (time(i,1)>5)&&(time(i,1) <9)
        time(i,4) = 2;
    elseif (time(i,1)>8)&&(time(i,1) <12)
        time(i,4) = 3;    
    end
end
clear i 

for i = 1:n_data1
    if (time(i,3) >=10)&&(time(i,3) <=16)
        time(i,5) = 1; % day : 1
    end
end
clear i 
for i = 1:n_data1
    if (time(i,3) >=20)||(time(i,3) <=5)
        time(i,5) = 2; % night : 2
    end
end
clear i 

%% 
result_day = zeros(4,10);
result_night = zeros(4,10);
result_inter = zeros(4,10);
result_total = zeros(4,10);
% high-medium-low urban and total, mean and standard error, 5 components
for i = 1:5 % components (from 3 to 7)
    % number
    n_day = zeros(4,1);
    n_night = zeros(4,1);
    n_inter = zeros(4,1);
    n_total = zeros(4,1);
    for j = 1:35
        for k = 1:n_data1 % 87649:n_data1
            if (data2((j-1)*n_data1+k,2+i)>0)&&(time(k,5)==1)
                n_day(4,1) = n_day(4,1)+1;
                n_day(i_stn(j,1),1) = n_day(i_stn(j,1),1)+1;
                n_total(4,1) = n_total(4,1)+1;
                n_total(i_stn(j,1),1) = n_total(i_stn(j,1),1)+1;
            elseif (data2((j-1)*n_data1+k,2+i)>0)&&(time(k,5)==2)
                n_night(4,1) = n_night(4,1)+1;
                n_night(i_stn(j,1),1) = n_night(i_stn(j,1),1)+1;
                n_total(4,1) = n_total(4,1)+1;
                n_total(i_stn(j,1),1) = n_total(i_stn(j,1),1)+1;
            elseif (data2((j-1)*n_data1+k,2+i)>0)&&(time(k,5)==0)
                n_inter(4,1) = n_inter(4,1)+1;
                n_inter(i_stn(j,1),1) = n_inter(i_stn(j,1),1)+1;
                n_total(4,1) = n_total(4,1)+1;
                n_total(i_stn(j,1),1) = n_total(i_stn(j,1),1)+1;
            end
        end
    end
    temp_day1 = zeros(n_day(1,1),1);
    temp_day2 = zeros(n_day(2,1),1);
    temp_day3 = zeros(n_day(3,1),1);
    temp_day4 = zeros(n_day(4,1),1);
    temp_night1 = zeros(n_night(1,1),1);
    temp_night2 = zeros(n_night(2,1),1);
    temp_night3 = zeros(n_night(3,1),1);
    temp_night4 = zeros(n_night(4,1),1);
    temp_inter1 = zeros(n_inter(1,1),1);
    temp_inter2 = zeros(n_inter(2,1),1);
    temp_inter3 = zeros(n_inter(3,1),1);
    temp_inter4 = zeros(n_inter(4,1),1);    
    temp_total1 = zeros(n_total(1,1),1);
    temp_total2 = zeros(n_total(2,1),1);
    temp_total3 = zeros(n_total(3,1),1);
    temp_total4 = zeros(n_total(4,1),1);
    n_day = zeros(4,1);
    n_night = zeros(4,1);
    n_inter = zeros(4,1);
    n_total = zeros(4,1);
    % sampling
    for j = 1:35
        for k = 1:n_data1 % 87649:n_data1
            if (data2((j-1)*n_data1+k,2+i)>0)&&(time(k,5)==1)
                n_day(4,1) = n_day(4,1)+1;
                n_day(i_stn(j,1),1) = n_day(i_stn(j,1),1)+1;
                n_total(4,1) = n_total(4,1)+1;
                n_total(i_stn(j,1),1) = n_total(i_stn(j,1),1)+1;
                if (i_stn(j,1)==1)
                    temp_day1(n_day(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_day4(n_day(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total1(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);
                elseif (i_stn(j,1)==2)
                    temp_day2(n_day(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_day4(n_day(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total2(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);                    
                elseif (i_stn(j,1)==3)
                    temp_day3(n_day(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_day4(n_day(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total3(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);                    
                end
            elseif (data2((j-1)*n_data1+k,2+i)>0)&&(time(k,5)==2)
                n_night(4,1) = n_night(4,1)+1;
                n_night(i_stn(j,1),1) = n_night(i_stn(j,1),1)+1;
                n_total(4,1) = n_total(4,1)+1;
                n_total(i_stn(j,1),1) = n_total(i_stn(j,1),1)+1;
                if (i_stn(j,1)==1)
                    temp_night1(n_night(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_night4(n_night(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total1(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);
                elseif (i_stn(j,1)==2)
                    temp_night2(n_night(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_night4(n_night(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total2(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);                    
                elseif (i_stn(j,1)==3)
                    temp_night3(n_night(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_night4(n_night(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total3(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);                    
                end                
            elseif (data2((j-1)*n_data1+k,2+i)>0)&&(time(k,5)==0)
                n_inter(4,1) = n_inter(4,1)+1;
                n_inter(i_stn(j,1),1) = n_inter(i_stn(j,1),1)+1;
                n_total(4,1) = n_total(4,1)+1;
                n_total(i_stn(j,1),1) = n_total(i_stn(j,1),1)+1;
                if (i_stn(j,1)==1)
                    temp_inter1(n_inter(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_inter4(n_inter(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total1(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);
                elseif (i_stn(j,1)==2)
                    temp_inter2(n_inter(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_inter4(n_inter(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total2(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);                    
                elseif (i_stn(j,1)==3)
                    temp_inter3(n_inter(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_inter4(n_inter(4,1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total3(n_total(i_stn(j,1),1),1) = data2((j-1)*n_data1+k,2+i);
                    temp_total4(n_total(4,1),1) = data2((j-1)*n_data1+k,2+i);                    
                end                
            end
        end
    end    
    % arrange results
    for j = 1:4
        if j==1
            result_day(j,(i-1)*2+1) = mean(temp_day1(:,1));
            result_night(j,(i-1)*2+1) = mean(temp_night1(:,1));
            result_inter(j,(i-1)*2+1) = mean(temp_inter1(:,1));
            result_total(j,(i-1)*2+1) = mean(temp_total1(:,1));
            result_day(j,(i-1)*2+2) = std(temp_day1(:,1))/(n_day(j,1))^(0.5);
            result_night(j,(i-1)*2+2) = std(temp_night1(:,1))/(n_night(j,1))^(0.5);
            result_inter(j,(i-1)*2+2) = std(temp_inter1(:,1))/(n_inter(j,1))^(0.5);
            result_total(j,(i-1)*2+2) = std(temp_total1(:,1))/(n_total(j,1))^(0.5);
        elseif j==2
            result_day(j,(i-1)*2+1) = mean(temp_day2(:,1));
            result_night(j,(i-1)*2+1) = mean(temp_night2(:,1));
            result_inter(j,(i-1)*2+1) = mean(temp_inter2(:,1));
            result_total(j,(i-1)*2+1) = mean(temp_total2(:,1));
            result_day(j,(i-1)*2+2) = std(temp_day2(:,1))/(n_day(j,1))^(0.5);
            result_night(j,(i-1)*2+2) = std(temp_night2(:,1))/(n_night(j,1))^(0.5);
            result_inter(j,(i-1)*2+2) = std(temp_inter2(:,1))/(n_inter(j,1))^(0.5);
            result_total(j,(i-1)*2+2) = std(temp_total2(:,1))/(n_total(j,1))^(0.5);
        elseif j==3
            result_day(j,(i-1)*2+1) = mean(temp_day3(:,1));
            result_night(j,(i-1)*2+1) = mean(temp_night3(:,1));
            result_inter(j,(i-1)*2+1) = mean(temp_inter3(:,1));
            result_total(j,(i-1)*2+1) = mean(temp_total3(:,1));
            result_day(j,(i-1)*2+2) = std(temp_day3(:,1))/(n_day(j,1))^(0.5);
            result_night(j,(i-1)*2+2) = std(temp_night3(:,1))/(n_night(j,1))^(0.5);
            result_inter(j,(i-1)*2+2) = std(temp_inter3(:,1))/(n_inter(j,1))^(0.5);
            result_total(j,(i-1)*2+2) = std(temp_total3(:,1))/(n_total(j,1))^(0.5);            
        elseif j==4
            result_day(j,(i-1)*2+1) = mean(temp_day4(:,1));
            result_night(j,(i-1)*2+1) = mean(temp_night4(:,1));
            result_inter(j,(i-1)*2+1) = mean(temp_inter4(:,1));
            result_total(j,(i-1)*2+1) = mean(temp_total4(:,1));
            result_day(j,(i-1)*2+2) = std(temp_day4(:,1))/(n_day(j,1))^(0.5);
            result_night(j,(i-1)*2+2) = std(temp_night4(:,1))/(n_night(j,1))^(0.5);
            result_inter(j,(i-1)*2+2) = std(temp_inter4(:,1))/(n_inter(j,1))^(0.5);
            result_total(j,(i-1)*2+2) = std(temp_total4(:,1))/(n_total(j,1))^(0.5);            
        end
    end
end
clear i j k

clear temp_day1 temp_day2 temp_day3 temp_day4
clear temp_night1 temp_night2 temp_night3 temp_night4
clear temp_inter1 temp_inter2 temp_inter3 temp_inter4
clear temp_total1 temp_total2 temp_total3 temp_total4
clear n_day n_night n_inter n_total
clear temp_n temp_sum