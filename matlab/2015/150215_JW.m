tic

clear all
clc

%% BOSUNG DATA PROCESS CODE
%   140710 Keunmin Lee; coding
%   140712 Keunmin Lee; gap checking process
%   150215 Je-Woo Hong; modify for new-version BS data

%% 10m
dataDir_10 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_10 = dir(fullfile(dataDir_10, '*_010M'))';
dataAll_10 = [];
for i=1:length(dataName_10)
    tempData_10 = importdata(fullfile(dataDir_10, dataName_10(i).name));
    dataAll_10 = [dataAll_10; tempData_10];
end
xlswrite('201501_010M.xls', dataAll_10)
clear dataDir_10 dataName_10 tempData_10 i

%% 20m
dataDir_20 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_20 = dir(fullfile(dataDir_20, '*_020M'))';
dataAll_20 = [];
for i=1:length(dataName_20)
    tempData_20 = importdata(fullfile(dataDir_20, dataName_20(i).name));
    dataAll_20 = [dataAll_20; tempData_20];
end
xlswrite('201501_020M.xls', dataAll_20)
clear dataDir_20 dataName_20 tempData_20 i

%% 40m
dataDir_40 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_40 = dir(fullfile(dataDir_40, '*_040M'))';
dataAll_40 = [];
for i=1:length(dataName_40)
    tempData_40 = importdata(fullfile(dataDir_40, dataName_40(i).name));
    dataAll_40 = [dataAll_40; tempData_40];
end
xlswrite('201501_040M.xls', dataAll_40)
clear dataDir_40 dataName_40 tempData_40 i

%% 60m
dataDir_60 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_60 = dir(fullfile(dataDir_60, '*_060M'))';
dataAll_60 = [];
for i=1:length(dataName_60)
    tempData_60 = importdata(fullfile(dataDir_60, dataName_60(i).name));
    dataAll_60 = [dataAll_60; tempData_60];
end
xlswrite('201501_060M.xls', dataAll_60)
clear dataDir_60 dataName_60 tempData_60 i

%% 80m
dataDir_80 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_80 = dir(fullfile(dataDir_80, '*_080M'))';
dataAll_80 = [];
for i=1:length(dataName_80)
    tempData_80 = importdata(fullfile(dataDir_80, dataName_80(i).name));
    dataAll_80 = [dataAll_80; tempData_80];
end
xlswrite('201501_080M.xls', dataAll_80)
clear dataDir_80 dataName_80 tempData_80 i

%% 100m
dataDir_100 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_100 = dir(fullfile(dataDir_100, '*_100M'))';
dataAll_100 = [];
for i=1:length(dataName_100)
    tempData_100 = importdata(fullfile(dataDir_100, dataName_100(i).name));
    dataAll_100 = [dataAll_100; tempData_100];
end
xlswrite('201501_100M.xls', dataAll_100)
clear dataDir_100 dataName_100 tempData_100 i

%% 140m
dataDir_140 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_140 = dir(fullfile(dataDir_140, '*_140M'))';
dataAll_140 = [];
for i=1:length(dataName_140)
    tempData_140 = importdata(fullfile(dataDir_140, dataName_140(i).name));
    dataAll_140 = [dataAll_140; tempData_140];
end
xlswrite('201501_140M.xls', dataAll_140)
clear dataDir_140 dataName_140 tempData_140 i

%% 180m
dataDir_180 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_180 = dir(fullfile(dataDir_180, '*_180M'))';
dataAll_180 = [];
for i=1:length(dataName_180)
    tempData_180 = importdata(fullfile(dataDir_180, dataName_180(i).name));
    dataAll_180 = [dataAll_180; tempData_180];
end
xlswrite('201501_180M.xls', dataAll_180)
clear dataDir_180 dataName_180 tempData_180 i

%% 220m
dataDir_220 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_220 = dir(fullfile(dataDir_220, '*_220M'))';
dataAll_220 = [];
for i=1:length(dataName_220)
    tempData_220 = importdata(fullfile(dataDir_220, dataName_220(i).name));
    dataAll_220 = [dataAll_220; tempData_220];
end
xlswrite('201501_220M.xls', dataAll_220)
clear dataDir_220 dataName_220 tempData_220 i

%% 260m
dataDir_260 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_260 = dir(fullfile(dataDir_260, '*_260M'))';
dataAll_260 = [];
for i=1:length(dataName_260)
    tempData_260 = importdata(fullfile(dataDir_260, dataName_260(i).name));
    dataAll_260 = [dataAll_260; tempData_260];
end
xlswrite('201501_260M.xls', dataAll_260)
clear dataDir_260 dataName_260 tempData_260 i

%% 300m
dataDir_300 = 'e:\2015_대학원\_보성타워\tower 1min\201501\'; % 폴더명
dataName_300 = dir(fullfile(dataDir_300, '*_300M'))';
dataAll_300 = [];
for i=1:length(dataName_300)
    tempData_300 = importdata(fullfile(dataDir_300, dataName_300(i).name));
    dataAll_300 = [dataAll_300; tempData_300];
end
xlswrite('201501_300M.xls', dataAll_300)
clear dataDir_300 dataName_300 tempData_300 i

%%
len = size(dataAll_10,1);
temperate = zeros(len,12);
humidity = zeros(len,12);
WS = zeros(len,12);
WD = zeros(len,12);
P = zeros(len,12);
clear len

for i=1:length(dataAll_10)
    temperate(i,1)=dataAll_10(i,3);
    humidity(i,1)=dataAll_10(i,4);
    WS(i,1)=dataAll_10(i,5);
    WD(i,1)=dataAll_10(i,6);
    P(i,1)=dataAll_10(i,9);
end
for i=1:length(dataAll_20)
    temperate(i,2)=dataAll_20(i,3);
    humidity(i,2)=dataAll_20(i,4);
    WS(i,2)=dataAll_20(i,5);
    WD(i,2)=dataAll_20(i,6);
    P(i,2)=dataAll_20(i,9);
end
for i=1:length(dataAll_40)
    temperate(i,3)=dataAll_40(i,3);
    humidity(i,3)=dataAll_40(i,4);
    WS(i,3)=dataAll_40(i,5);
    WD(i,3)=dataAll_40(i,6);
    P(i,3)=dataAll_40(i,9);
end
for i=1:length(dataAll_60)
    temperate(i,4)=dataAll_60(i,3);
    humidity(i,4)=dataAll_60(i,4);
    WS(i,4)=dataAll_60(i,5);
    WD(i,4)=dataAll_60(i,6);
    P(i,4)=dataAll_60(i,9);
end
for i=1:length(dataAll_80)
    temperate(i,5)=dataAll_80(i,3);
    humidity(i,5)=dataAll_80(i,4);
    WS(i,5)=dataAll_80(i,5);
    WD(i,5)=dataAll_80(i,6);
    P(i,5)=dataAll_80(i,9);
end
for i=1:length(dataAll_100)
    temperate(i,6)=dataAll_100(i,3);
    humidity(i,6)=dataAll_100(i,4);
    WS(i,6)=dataAll_100(i,5);
    WD(i,6)=dataAll_100(i,6);
    P(i,6)=dataAll_100(i,9);
end
for i=1:length(dataAll_140)
    temperate(i,7)=dataAll_140(i,3);
    humidity(i,7)=dataAll_140(i,4);
    WS(i,7)=dataAll_140(i,5);
    WD(i,7)=dataAll_140(i,6);
    P(i,7)=dataAll_140(i,9);
end
for i=1:length(dataAll_180)
    temperate(i,8)=dataAll_180(i,3);
    humidity(i,8)=dataAll_180(i,4);
    WS(i,8)=dataAll_180(i,5);
    WD(i,8)=dataAll_180(i,6);
    P(i,8)=dataAll_180(i,9);
end
for i=1:length(dataAll_220)
    temperate(i,9)=dataAll_220(i,3);
    humidity(i,9)=dataAll_220(i,4);
    WS(i,9)=dataAll_220(i,5);
    WD(i,9)=dataAll_220(i,6);
    P(i,9)=dataAll_220(i,9);
end
for i=1:length(dataAll_260)
    temperate(i,10)=dataAll_260(i,3);
    humidity(i,10)=dataAll_260(i,4);
    WS(i,10)=dataAll_260(i,5);
    WD(i,10)=dataAll_260(i,6);
    P(i,10)=dataAll_260(i,9);
end
for i=1:length(dataAll_300)
    temperate(i,11)=dataAll_300(i,3);
    humidity(i,11)=dataAll_300(i,4);
    WS(i,11)=dataAll_300(i,5);
    WD(i,11)=dataAll_300(i,6);
    P(i,11)=dataAll_300(i,9);
end
for i=1:length(dataAll_300)
    temperate(i,12)=dataAll_300(i,2);
    humidity(i,12)=dataAll_300(i,2);
    WS(i,12)=dataAll_300(i,2);
    WD(i,12)=dataAll_300(i,2);
    P(i,12)=dataAll_300(i,2);
end
clear i 

xlswrite('201501_T.xls', temperate)
xlswrite('201501_RH.xls', humidity)
xlswrite('201501_WD.xls', WD)
xlswrite('201501_WS.xls', WS)
xlswrite('201501_P.xls', P)

toc

