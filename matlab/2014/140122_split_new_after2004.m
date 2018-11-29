%% before run
% 1. import "raw_data" as "raw"
% 2. change "size" in this code
% 3. change year (ex 2004-> your yr)


%% split 2004 (new version)
% size = 3711260; % 03
size = 3747929; % 04
% size = 3760382 % 05
% size = 3848414 % 09
% size = 3879104 % 10
% size = 3933838 %11
% size = 3977019 %12

data = raw;
% wait = waitbar(0,'...split...');
for i = 1:size
%     waitbar((i/(size)),wait,sprintf('%f',i/(size)*100))

%     data(i,1) = strrep(data(i, 1),'"	"',',');
%     data(i,1) = strrep(data(i, 1),'"',',');
%     data(i,1) = strrep(data(i,1),',,',',-999,');
%     data(i,1) = strrep(data(i,1),',,',',');
    data(i,1) = strrep(data(i, 1),'""','"-999"');
end
i1 = i;
i1
clear i i1
data = regexp(data,',','split');

aws2004 = zeros(size-1,12);
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
for i = 2:size
   for j = 5:8
        temp1 = data{i,1}(1,j);
        aws2004(i-1,j+4) = str2num(temp1{1,1});
   end
   temp1 = regexp(data{i,1}{1,1},'	','split'); 
   aws2004(i-1,1) = str2num(temp1{1,1});
   aws2004(i-1,2) = str2num(temp1{1,2});
   
   temp1 = regexp(data{i,1}{1,3},'	','split'); 
   aws2004(i-1,3) = str2num(temp1{1,2});
   aws2004(i-1,4) = str2num(temp1{1,3});
   
   temp1 = regexp(data{i,1}{1,4},' ','split');
   aws2004(i-1,8) = str2num(temp1{1,2});
   temp2 = regexp(temp1{1,1},'-','split');
   aws2004(i-1,5) = str2num(temp2{1,1});
   aws2004(i-1,6) = str2num(temp2{1,2});
   aws2004(i-1,7) = str2num(temp2{1,3});
end
i2 = i;
i2
clear i i2 j temp1 temp2


%% Station number
stn_n = 0;
stn_id = 0;
temp = 0;
for i = 1:size-1
    if temp ~= aws2004(i,2)
        stn_n = stn_n + 1;
        temp = aws2004(i,2);
        stn_id(stn_n,1) = stn_n;
        % order
        stn_id(stn_n,2) = temp;
        % stn_id
    end
end
clear i temp


%% GEV analysis
GEV = zeros(stn_n,4);
% 1 : order 
% 2 : stn_id
% 3 : daily max rainfall
% 4 : 5-days max rainfall
temp = 0; % for stn_id
temp1 = 0; % for order
% for daily max rainfall
temp2_2 = 0; % for DOY
temp2_3 = 0; % for calculation
temp2_4 = 0; % for calculation
% for 5-days max rainfall
temp3_2 = zeros(5,1); % for number of days
temp3_3 = 0; % for calculation

for i = 1:size-1
    if temp ~= aws2004(i,2) % changes in stn_id
        temp = aws2004(i,2);
        temp1 = temp1 + 1;
        temp2_2 = aws2004(i,7);
        GEV(temp1,1) = temp1;
        GEV(temp1,2) = temp;
        if temp1 >= 2
            GEV(temp1-1,3) = temp2_4;
            temp2_4 = 0;
            GEV(temp1-1,4) = temp3_3;
            temp3_3 = 0;
        end
    end
    if temp2_2 ~= aws2004(i,7) % changes in day within same stn_id
        temp2_2 = aws2004(i,7);
        if temp2_3 > temp2_4
            temp2_4 = temp2_3;
        end
        temp2_3 = aws2004(i,10);

        if temp3_2(5,1) > temp3_3
            temp3_3 = temp3_2(5,1);
        end
        for j = 5:-1:2
            temp3_2(j,1) = temp3_2(j-1,1);
        end
        temp3_2(1,1) = 0;
    end
    if aws2004(i,10)>0 && aws2004(i,10)<3000
        temp2_3 = temp2_3 + aws2004(i,10);
        for j = 1:5
            temp3_2(j,1) = temp3_2(j,1) + aws2004(i,10);
        end
    end
    if i == size-1 % for the last stn
        GEV(temp1,3) = temp2_4;
        GEV(temp1,4) = temp3_3;
    end
end
clear i j temp temp1 temp2_2 temp2_3 temp2_4 temp3_2 temp3_3


%% extract GEV for Seoul
aws_n_seoul = [108,7;110,20;116,15;400,4;401,4;402,1;403,1;404,2;405,10;406,16;407,12;408,11;409,3;410,9;411,8;412,18;413,5;414,12;415,1;416,18;417,3;418,21;419,15;421,19;422,18;423,6;424,4;425,18;437,4;453,14;492,13;509,17;510,4;590,20;];
stn_n_seoul = 34;
temp = 0; % for address
GEV_seoul = zeros(stn_n_seoul,5);
% 1 : order 
% 2 : stn_id
% 3 : LCZ
% 4 : daily max rainfall
% 5 : 5-days max rainfall
for i = 1:stn_n_seoul
    GEV_seoul(i,1) = i;
    GEV_seoul(i,2) = aws_n_seoul(i,1);
    GEV_seoul(i,3) = aws_n_seoul(i,2);
    temp = find(GEV_seoul(i,2)==GEV(:,2));
    if temp>0
        GEV_seoul(i,4) = GEV(temp,3);
        GEV_seoul(i,5) = GEV(temp,4);
    end
    temp = 0;
end
clear i temp 


%% finish
clear aws_n_seoul stn_id stn_n stn_n_seoul size data
% close(wait);
clear wait


