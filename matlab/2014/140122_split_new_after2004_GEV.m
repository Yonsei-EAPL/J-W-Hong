%% before run
% 1. import "raw_data" as "raw"
% 2. change "n_size" in this code
% 3. change year (ex 2012-> your yr)


%% split 2012 (new version)
% n_size = 3711260 % 03
% n_size = 3747929 % 04
% n_size = 3760382 % 05
% n_size = 3848414 % 09
% n_size = 3879104 % 10
% n_size = 3933838 %11
n_size = 3977019 %12

data = raw;
% wait = waitbar(0,'...split...');
for i = 1:n_size
%     waitbar((i/(n_size)),wait,sprintf('%f',i/(n_size)*100))

%     data(i,1) = strrep(data(i, 1),'"	"',',');
%     data(i,1) = strrep(data(i, 1),'"',',');
%     data(i,1) = strrep(data(i,1),',,',',-999,');
%     data(i,1) = strrep(data(i,1),',,',',');
    data(i,1) = strrep(data(i, 1),'""','"-999"');
    data(i,1) = strrep(data(i, 1),' ','');
    data(i,1) = strrep(data(i, 1),'"','');
end
i1 = i;
i1
clear i i1
data = regexp(data,',','split');

aws2012 = zeros(n_size-1,12);
% 1 : num
% 2 : stn_num
% 3 : lat
% 4 : lon
% 5 : year
% 6 : month
% 7 : day
% 8 : hour
% 9 : Temp
% 10 : Rainfall
% 11 : Wind-speed
% 12 : Wind-direction
for i = 2:n_size
   for j = 7:10
       aws2012(i-1,j+2) = str2num(data{i,1}{1,j});
   end
%    temp1 = regexp(data{i,1}{1,1},'	','split'); 
%    aws2012(i-1,1) = str2num(temp1{1,1});
%    aws2012(i-1,2) = str2num(temp1{1,2});
   aws2012(i-1,1) = str2num(data{i,1}{1,1});
   aws2012(i-1,2) = str2num(data{i,1}{1,2});
   
%    temp1 = regexp(data{i,1}{1,3},'	','split'); 
%    aws2012(i-1,3) = str2num(temp1{1,2});
%    aws2012(i-1,4) = str2num(temp1{1,3});
   aws2012(i-1,3) = str2num(data{i,1}{1,4});
   aws2012(i-1,4) = str2num(data{i,1}{1,5});   

   temp1 = regexp(data{i,1}{1,6},' ','split');
   temp2 = regexp(temp1{1,1},'-','split');
   aws2012(i-1,5) = str2num(temp2{1,1});
   aws2012(i-1,6) = str2num(temp2{1,2});
   aws2012(i-1,7) = str2num(temp2{1,3});
   aws2012(i-1,8) = mod(aws2012(i-1,7),100);
   temp = mod(aws2012(i-1,7),100);
   aws2012(i-1,7) = (aws2012(i-1,7) - temp)/100;
end
i2 = i;
i2
clear i i2 j temp1 temp2


%% Station number
stn_n = 0;
stn_id = 0;
temp = 0;
for i = 1:n_size-1
    if temp ~= aws2012(i,2)
        stn_n = stn_n + 1;
        temp = aws2012(i,2);
        stn_id(stn_n,1) = stn_n;
        % order
        stn_id(stn_n,2) = temp;
        % stn_id
    end
end
clear i temp


%% GEV analysis
GEV = zeros(stn_n,3);
% 1 : order 
% 2 : stn_id
% 3 : max rainfall
for i = 1:stn_n
    GEV(i,1) = i;
    GEV(i,2) = stn_id(i,2);
    temp = find(aws2012(:,2)==stn_id(i,2));
    temp2 = size(temp);
    temp1 = zeros(temp2(1,1),1);
    for j = 1:temp2(1,1)
        temp1(j,1) = aws2012(temp(j,1),10);
    end
    temp3 = max(temp1);
    GEV(i,3) = temp3;
end
clear i j temp temp1 temp2 temp3


%% extract GEV for Seoul
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
    temp1 = find(aws2012(:,2)==aws_n_seoul(i,1));
    temp2 = size(temp1);
    if temp2(1,1)>0
        for j = 1:temp2(1,1)
            if aws2012(temp1(j,1),6) > 1
                temp3 =0;
                for k = 1:(aws2012(temp1(j,1),6)-1)
                    temp3 = temp3 + n_month(k,1);
                end
                temp3 = temp3*24 + (aws2012(temp1(j,1),7)-1)*24 + aws2012(temp1(j,1),8)+1;
                temp(temp3,i+1) = aws2012(temp1(j,1),9);
            else
                temp3 = (aws2012(temp1(j,1),7)-1)*24 + aws2012(temp1(j,1),8) +1;
                temp(temp3,i+1) = aws2012(temp1(j,1),9);
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


