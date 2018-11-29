parpool(2);
dataDir = 'H:\b1_EAPL\b1_JW_Observation\jj2018\csv' %'E:\EAPL\JW_Observation\은평뉴타운\raw_csv'; % folder name % EP_1day_csv
dataName = dir(fullfile(dataDir, 'CSV_8201.ts_data_*'))'; % 
n= length(dataName);
data = cell(n,1);
parfor i = 13857:n
    data{i,1} = importdata(fullfile(dataDir, dataName(i).name));
    dataName(i).name = strcat('matlab_',dataName(i).name);
    data{i,1}(:,14) = data{i,1}(:,14) - data{i,1}(:,15);
    csvwrite(dataName(i).name, data{i,1});
    data{i,1}=0;
end
delete(gcp('nocreate'))