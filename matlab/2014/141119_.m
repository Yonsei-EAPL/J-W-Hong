%% BOSUNG DATA PROCESS CODE
%   140710 Keunmin Lee; coding
%   140712 Keunmin Lee; gap checking process

%% 10m
dataDir_10 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\10M\'; % 폴더명
dataName_10 = dir(fullfile(dataDir_10, '01MN-*'))';
dataAll_10 = [];
for i=1:length(dataName_10)
    tempData_10 = importdata(fullfile(dataDir_10, dataName_10(i).name));
    dataAll_10 = [dataAll_10; tempData_10];
end
clear dataDir_10 dataName_10 tempData_10 i

for i =1:length(dataAll_10)-1
    interval(i,1)=dataAll_10(i+1,1)-dataAll_10(i,1);
end
clear i
for  i=1:length(dataAll_10)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_10,2)
                gap(x,y)=999;
            end
        end
        dataAll_10=[dataAll_10(1:i,:);gap(1:end,:); dataAll_10(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_10);
    dataAll_10(i,2)=(dataAll_10(i,2)-10000)/100;
    dataAll_10(i,5)=dataAll_10(i,5)/10;
    dataAll_10(i,8)=dataAll_10(i,8)/10;
    dataAll_10(i,10)=round(dataAll_10(i,10)/10);
end

clear i 

%% 20m
dataDir_20 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\20M\'; % 폴더명
dataName_20 = dir(fullfile(dataDir_20, '01MN-*'))';
dataAll_20 = [];
for i=1:length(dataName_20)
    tempData_20 = importdata(fullfile(dataDir_20, dataName_20(i).name));
    dataAll_20 = [dataAll_20; tempData_20];
end
clear dataDir_20 dataName_20 tempData_20 i

for i =1:length(dataAll_20)-1
    interval(i,1)=dataAll_20(i+1,1)-dataAll_20(i,1);
end
clear i
for  i=1:length(dataAll_20)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_20,2)
                gap(x,y)=999;
            end
        end
        dataAll_20=[dataAll_20(1:i,:);gap(1:end,:); dataAll_20(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_20);
    dataAll_20(i,2)=(dataAll_20(i,2)-10000)/100;
    dataAll_20(i,5)=dataAll_20(i,5)/10;
    dataAll_20(i,8)=dataAll_20(i,8)/10;
    dataAll_20(i,10)=round(dataAll_20(i,10)/10);
end

clear i 

%% 40m
dataDir_40 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\40M\'; % 폴더명
dataName_40 = dir(fullfile(dataDir_40, '01MN-*'))';
dataAll_40 = [];
for i=1:length(dataName_40)
    tempData_40 = importdata(fullfile(dataDir_40, dataName_40(i).name));
    dataAll_40 = [dataAll_40; tempData_40];
end
clear dataDir_40 dataName_40 tempData_40 i

for i =1:length(dataAll_40)-1
    interval(i,1)=dataAll_40(i+1,1)-dataAll_40(i,1);
end
clear i
for  i=1:length(dataAll_40)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_40,2)
                gap(x,y)=999;
            end
        end
        dataAll_40=[dataAll_40(1:i,:);gap(1:end,:); dataAll_40(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_40);
    dataAll_40(i,2)=(dataAll_40(i,2)-10000)/100;
    dataAll_40(i,5)=dataAll_40(i,5)/10;
    dataAll_40(i,8)=dataAll_40(i,8)/10;
    dataAll_40(i,10)=round(dataAll_40(i,10)/10);
end

clear i 

%% 60m
dataDir_60 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\60M\'; % 폴더명
dataName_60 = dir(fullfile(dataDir_60, '01MN-*'))';
dataAll_60 = [];
for i=1:length(dataName_60)
    tempData_60 = importdata(fullfile(dataDir_60, dataName_60(i).name));
    dataAll_60 = [dataAll_60; tempData_60];
end
clear dataDir_60 dataName_60 tempData_60 i

for i =1:length(dataAll_60)-1
    interval(i,1)=dataAll_60(i+1,1)-dataAll_60(i,1);
end
clear i
for  i=1:length(dataAll_60)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_60,2)
                gap(x,y)=999;
            end
        end
        dataAll_60=[dataAll_60(1:i,:);gap(1:end,:); dataAll_60(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_60);
    dataAll_60(i,2)=(dataAll_60(i,2)-10000)/100;
    dataAll_60(i,5)=dataAll_60(i,5)/10;
    dataAll_60(i,8)=dataAll_60(i,8)/10;
    dataAll_60(i,10)=round(dataAll_60(i,10)/10);
end
clear i 

%% 80m
dataDir_80 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\80M\'; % 폴더명
dataName_80 = dir(fullfile(dataDir_80, '01MN-*'))';
dataAll_80 = [];
for i=1:length(dataName_80)
    tempData_80 = importdata(fullfile(dataDir_80, dataName_80(i).name));
    dataAll_80 = [dataAll_80; tempData_80];
end
clear dataDir_80 dataName_80 tempData_80 i

for i =1:length(dataAll_80)-1
    interval(i,1)=dataAll_80(i+1,1)-dataAll_80(i,1);
end
clear i
for  i=1:length(dataAll_80)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_80,2)
                gap(x,y)=999;
            end
        end
        dataAll_80=[dataAll_80(1:i,:);gap(1:end,:); dataAll_80(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_80);
    dataAll_80(i,2)=(dataAll_80(i,2)-10000)/100;
    dataAll_80(i,5)=dataAll_80(i,5)/10;
    dataAll_80(i,8)=dataAll_80(i,8)/10;
    dataAll_80(i,10)=round(dataAll_80(i,10)/10);
end
clear i 

%% 100m
dataDir_100 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\100M\'; % 폴더명
dataName_100 = dir(fullfile(dataDir_100, '01MN-*'))';
dataAll_100 = [];
for i=1:length(dataName_100)
    tempData_100 = importdata(fullfile(dataDir_100, dataName_100(i).name));
    dataAll_100 = [dataAll_100; tempData_100];
end
clear dataDir_100 dataName_100 tempData_100 i

for i =1:length(dataAll_100)-1
    interval(i,1)=dataAll_100(i+1,1)-dataAll_100(i,1);
end
clear i
for  i=1:length(dataAll_100)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_100,2)
                gap(x,y)=999;
            end
        end
        dataAll_100=[dataAll_100(1:i,:);gap(1:end,:); dataAll_100(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_100);
    dataAll_100(i,2)=(dataAll_100(i,2)-10000)/100;
    dataAll_100(i,5)=dataAll_100(i,5)/10;
    dataAll_100(i,8)=dataAll_100(i,8)/10;
    dataAll_100(i,10)=round(dataAll_100(i,10)/10);
end
clear i 

%% 140m
dataDir_140 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\140M\'; % 폴더명
dataName_140 = dir(fullfile(dataDir_140, '01MN-*'))';
dataAll_140 = [];
for i=1:length(dataName_140)
    tempData_140 = importdata(fullfile(dataDir_140, dataName_140(i).name));
    dataAll_140 = [dataAll_140; tempData_140];
end
clear dataDir_140 dataName_140 tempData_140 i

for i =1:length(dataAll_140)-1
    interval(i,1)=dataAll_140(i+1,1)-dataAll_140(i,1);
end
clear i
for  i=1:length(dataAll_140)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_140,2)
                gap(x,y)=999;
            end
        end
        dataAll_140=[dataAll_140(1:i,:);gap(1:end,:); dataAll_140(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_140);
    dataAll_140(i,2)=(dataAll_140(i,2)-10000)/100;
    dataAll_140(i,5)=dataAll_140(i,5)/10;
    dataAll_140(i,8)=dataAll_140(i,8)/10;
    dataAll_140(i,10)=round(dataAll_140(i,10)/10);
end
clear i 

%% 180m
dataDir_180 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\180M\'; % 폴더명
dataName_180 = dir(fullfile(dataDir_180, '01MN-*'))';
dataAll_180 = [];
for i=1:length(dataName_180)
    tempData_180 = importdata(fullfile(dataDir_180, dataName_180(i).name));
    dataAll_180 = [dataAll_180; tempData_180];
end
clear dataDir_180 dataName_180 tempData_180 i

for i =1:length(dataAll_180)-1
    interval(i,1)=dataAll_180(i+1,1)-dataAll_180(i,1);
end
clear i
for  i=1:length(dataAll_180)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_180,2)
                gap(x,y)=999;
            end
        end
        dataAll_180=[dataAll_180(1:i,:);gap(1:end,:); dataAll_180(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_180);
    dataAll_180(i,2)=(dataAll_180(i,2)-10000)/100;
    dataAll_180(i,5)=dataAll_180(i,5)/10;
    dataAll_180(i,8)=dataAll_180(i,8)/10;
    dataAll_180(i,10)=round(dataAll_180(i,10)/10);
end
clear i 

%% 220m
dataDir_220 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\220M\'; % 폴더명
dataName_220 = dir(fullfile(dataDir_220, '01MN-*'))';
dataAll_220 = [];
for i=1:length(dataName_220)
    tempData_220 = importdata(fullfile(dataDir_220, dataName_220(i).name));
    dataAll_220 = [dataAll_220; tempData_220];
end
clear dataDir_220 dataName_220 tempData_220 i

for i =1:length(dataAll_220)-1
    interval(i,1)=dataAll_220(i+1,1)-dataAll_220(i,1);
end
clear i
for  i=1:length(dataAll_220)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_220,2)
                gap(x,y)=999;
            end
        end
        dataAll_220=[dataAll_220(1:i,:);gap(1:end,:); dataAll_220(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_220);
    dataAll_220(i,2)=(dataAll_220(i,2)-10000)/100;
    dataAll_220(i,5)=dataAll_220(i,5)/10;
    dataAll_220(i,8)=dataAll_220(i,8)/10;
    dataAll_220(i,10)=round(dataAll_220(i,10)/10);
end
clear i 

%% 260m
dataDir_260 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\260M\'; % 폴더명
dataName_260 = dir(fullfile(dataDir_260, '01MN-*'))';
dataAll_260 = [];
for i=1:length(dataName_260)
    tempData_260 = importdata(fullfile(dataDir_260, dataName_260(i).name));
    dataAll_260 = [dataAll_260; tempData_260];
end
clear dataDir_260 dataName_260 tempData_260 i

for i =1:length(dataAll_260)-1
    interval(i,1)=dataAll_260(i+1,1)-dataAll_260(i,1);
end
clear i
for  i=1:length(dataAll_260)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_260,2)
                gap(x,y)=999;
            end
        end
        dataAll_260=[dataAll_260(1:i,:);gap(1:end,:); dataAll_260(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_260);
    dataAll_260(i,2)=(dataAll_260(i,2)-10000)/100;
    dataAll_260(i,5)=dataAll_260(i,5)/10;
    dataAll_260(i,8)=dataAll_260(i,8)/10;
    dataAll_260(i,10)=round(dataAll_260(i,10)/10);
    dataAll_260(i,12)=dataAll_260(i,12)/100;
    dataAll_260(i,13)=round(dataAll_260(i,13)/10);
    dataAll_260(i,15)=round(dataAll_260(i,15)/10);
end
clear i 

%% 300m
dataDir_300 = 'g:\외장_MATLAB\2014_대학원\141119_contour_보성\1분_기상연 자료\TWR_DATA_20141008_20141109\300M\'; % 폴더명
dataName_300 = dir(fullfile(dataDir_300, '01MN-*'))';
dataAll_300 = [];
for i=1:length(dataName_300)
    tempData_300 = importdata(fullfile(dataDir_300, dataName_300(i).name));
    dataAll_300 = [dataAll_300; tempData_300];
end
clear dataDir_300 dataName_300 tempData_300 i

for i =1:length(dataAll_300)-1
    interval(i,1)=dataAll_300(i+1,1)-dataAll_300(i,1);
end
clear i
for  i=1:length(dataAll_300)-1;
    if ((interval(i,1)>1)&&(interval(i,1)~=41)&&(interval(i,1)~=7641)&&(interval(i,1)~=697641)&&(interval(i,1)~=707641)&&(interval(i,1)~=88697641))
        for x=1:(interval(i,1)-1)
            for y=1:size(dataAll_300,2)
                gap(x,y)=999;
            end
        end
        dataAll_300=[dataAll_300(1:i,:);gap(1:end,:); dataAll_300(i+1:end,:)];
        clear x y gap
      end
end
clear i interval 

for i=1:length(dataAll_300);
    dataAll_300(i,2)=(dataAll_300(i,2)-10000)/100;
    dataAll_300(i,5)=dataAll_300(i,5)/10;
    dataAll_300(i,8)=dataAll_300(i,8)/10;
    dataAll_300(i,10)=round(dataAll_300(i,10)/10);
end
clear i 
toc
tic
%%