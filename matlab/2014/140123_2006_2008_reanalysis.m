clear result

%% using aws2008(n,11), the former data
% n_size = 5373790; % 06
% n_size = 5373790; % 07
n_size = 5596499; % 08
% n_size = ; % 12

aws2008_n = zeros(n_size,8);
for i = 1:n_size
    for j = 1:8
        aws2008_n(i,j) = -999;
    end
end
clear i j 
for i = 1:n_size
    aws2008_n(i,1) = i;
    aws2008_n(i,2) = aws2008(i,2);
    
    aws2008_n(i,3) = 2008;

    temp_m = mod(aws2008(i,1),1000000) - mod(aws2008(i,1),10000);
    aws2008_n(i,4) = temp_m/10000;
    temp_d = mod(aws2008(i,1),10000) - mod(aws2008(i,1),100);
    aws2008_n(i,5) = temp_d/100;
    temp_h = mod(aws2008(i,1),100);
    aws2008_n(i,6) = temp_h;
    
    aws2008_n(i,7) = aws2008(i,3);
    aws2008_n(i,8) = aws2008(i,4);    
end
clear i temp_m temp_d temp_h

aws2008 = aws2008_n;
clear aws2008_n


%% GEV analysis 

% Station number
stn_n = 0;
stn_id = 0;
temp = 0;
a = sort(aws2008(:,2));
for i = 1:n_size
    if temp ~= a(i,1)
        stn_n = stn_n + 1;
        temp = a(i,1);
        stn_id(stn_n,1) = stn_n;
        % order
        stn_id(stn_n,2) = temp;
        % stn_id
    end
end
clear i temp a

% GEV analysis
GEV = zeros(stn_n,3);
    % 1 : order 
    % 2 : stn_id
    % 3 : max rainfall
for i = 1:stn_n
    GEV(i,1) = i;
    GEV(i,2) = stn_id(i,2);
    temp = find(aws2008(:,2)==stn_id(i,2));
    temp2 = size(temp);
    temp1 = zeros(temp2(1,1),1);
    for j = 1:temp2(1,1)
        temp1(j,1) = aws2008(temp(j,1),8);
    end
    temp3 = max(temp1);
    GEV(i,3) = temp3;
end
clear i j temp temp1 temp2 temp3

% extract GEV for Seoul
aws_n_seoul = [108,7;110,20;116,15;400,4;401,4;402,1;403,1;404,2;405,10;406,16;407,12;408,11;409,3;410,9;411,8;412,18;413,5;414,12;415,1;416,18;417,3;418,21;419,15;421,19;422,18;423,6;424,4;425,18;437,4;453,14;492,13;509,17;510,4;590,20;];
stn_n_seoul = 34;
temp = 0; % for address
GEV_seoul = zeros(stn_n_seoul,5);
% 1 : order 
% 2 : stn_id
% 3 : LCZ
% 4 : max rainfall
% 5 : max UHI
for i = 1:stn_n_seoul
    GEV_seoul(i,1) = i;
    GEV_seoul(i,2) = aws_n_seoul(i,1);
    GEV_seoul(i,3) = aws_n_seoul(i,2);
    temp = find(GEV_seoul(i,2)==GEV(:,2));
    if temp>0
        GEV_seoul(i,4) = GEV(temp,3);
    end
    temp = 0;
end
clear i temp 

n_month = [31;29;31;30;31;30;31;31;30;31;30;31]; % for '04
% n_month = [31;28;31;30;31;30;31;31;30;31;30;31]; % for other yr
temp = zeros(sum(n_month)*24,stn_n_seoul+1); % for UHI
for i = 1:sum(n_month)*24
    temp(i,1) = i;
    for j = 1:stn_n_seoul
        temp(i,j+1) = -999;
    end
end
clear i j
for i = 1:stn_n_seoul
    temp1 = find(aws2008(:,2)==aws_n_seoul(i,1));
    temp2 = size(temp1);
    if temp2(1,1)>0
        for j = 1:temp2(1,1)
            if aws2008(temp1(j,1),4) > 1
                temp3 =0;
                for k = 1:(aws2008(temp1(j,1),4)-1)
                    temp3 = temp3 + n_month(k,1);
                end
                temp3 = temp3*24 + (aws2008(temp1(j,1),5)-1)*24 + aws2008(temp1(j,1),6)+1;
                temp(temp3,i+1) = aws2008(temp1(j,1),7);
            else
                temp3 = (aws2008(temp1(j,1),5)-1)*24 + aws2008(temp1(j,1),6) +1;
                temp(temp3,i+1) = aws2008(temp1(j,1),7);
            end
        end
    end
end
clear i j k  temp1 temp2 temp3

n_ref = 32; % for Gwan-ak site (509)
for i = 1:sum(n_month)*24
    if (temp(i,n_ref+1)>-20) && (temp(i,n_ref+1)<40)
        for j = 1:stn_n_seoul
            if (temp(i,j+1)>-20) && (temp(i,j+1)<40)
                temp(i,j+1) = temp(i,j+1) - temp(i,n_ref+1);
            else
                temp(i,j+1) = -999;
            end
        end
    else
        for j = 1:stn_n_seoul
            temp(i,j+1) = -999;
        end
    end
end
clear i j

for i = 1:stn_n_seoul
    GEV_seoul(i,5) = max(temp(:,i+1));
end
clear i

%% finish
clear aws_n_seoul stn_id stn_n stn_n_seoul n_size data n_ref rn n_month temp
% close(wait);
clear wait


