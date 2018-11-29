dataDir = 'H:\b1_EAPL\b1_JW_Observation\ep2018\_csv1' %'E:\EAPL\JW_Observation\은평뉴타운\raw_csv'; % folder name % EP_1day_csv
dataName = dir(fullfile(dataDir, 'CSV_3302.ts_data_*'))'; % 

for t=1:5000%length(dataName)
    data = importdata(fullfile(dataDir, dataName(t).name));  
    data(1,20)=0;
    filename = strcat('m',dataName(t).name);
    data(:,20) = data(:,14)-data(:,17);
    
    dlmwrite(filename, data); 
    clear data
end
clear t

