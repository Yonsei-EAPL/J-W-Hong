dataDir = 'E:\2016_대학원\_보성타워_16\1min\twr\2016\10\'; % folder
dataNameT = dir(fullfile(dataDir,'*_T*'));
dataNameRH = dir(fullfile(dataDir,'*_H*'));
dataNameWS = dir(fullfile(dataDir,'*_V*'));
dataNameWD = dir(fullfile(dataDir,'*_WD*'));

dataAllT = [];
for i=1:length(dataNameT)
    tempDataT = importdata(fullfile(dataDir, dataNameT(i).name));
    dataAllT = [dataAllT; tempDataT];
end
xlswrite('201610T.xls', dataAllT)

dataAllRH = [];
for i=1:length(dataNameRH)
    tempDataRH = importdata(fullfile(dataDir, dataNameRH(i).name));
    dataAllRH = [dataAllRH; tempDataRH];
end
xlswrite('201610RH.xls', dataAllRH)

dataAllWS = [];
for i=1:length(dataNameWS)
    tempDataWS = importdata(fullfile(dataDir, dataNameWS(i).name));
    dataAllWS = [dataAllWS; tempDataWS];
end
xlswrite('201610WS.xls', dataAllWS)

dataAllWD = [];
for i=1:length(dataNameWD)
    tempDataWD = importdata(fullfile(dataDir, dataNameWD(i).name));
    dataAllWD = [dataAllWD; tempDataWD];
end
xlswrite('201610WD.xls', dataAllWD)
