dataDir = 'H:\b1_EAPL\b1_JW_Observation\jj2018\csv' %'E:\EAPL\JW_Observation\은평뉴타운\raw_csv'; % folder name % EP_1day_csv
dataName = dir(fullfile(dataDir, 'CSV_8201.ts_data_*'))'; % 

for t=13621:length(dataName)
    data = importdata(fullfile(dataDir, dataName(t).name));  
    filename = strcat('matlab_',dataName(t).name);
    data(:,14) = data(:,14)-data(:,15);
    
    csvwrite(filename, data);
    clear data
end
clear t

