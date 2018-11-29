dataDir = 'E:\다운로드\mat_files (2)'; 
dataName = dir(fullfile(dataDir, '*.mat'))'; 

for t=1:length(dataName) 

    Data = load(fullfile(dataDir, dataName(t).name));
    DataField = fieldnames(Data);
    name = erase(dataName(t).name,".mat");
    name = [name '.txt'];
    dlmwrite(name, Data.(DataField{1}));
      
    
end
clear t

