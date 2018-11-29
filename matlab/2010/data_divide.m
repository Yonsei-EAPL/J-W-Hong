data_set=importdata(input('filename= ','s'));
site=input('site= ','s');
year=input('year= ','s');

data_size=size(data_set.data,1);

if data_size == 17760
    for i=1:12
        filename=strcat(site,'_',year,'_',num2str(i),'.csv');
    if i==1
        csvwrite(filename,data_set.data(1:1680,:));
    elseif i==2
        csvwrite(filename,data_set.data(1489:3072,:));
    elseif i==3 
        csvwrite(filename,data_set.data(2881:4560,:));
    elseif i==4
        csvwrite(filename,data_set.data(4369:6000,:));
    elseif i==5
        csvwrite(filename,data_set.data(5809:7488,:));
    elseif i==6
        csvwrite(filename,data_set.data(7297:8928,:));
    elseif i==7
        csvwrite(filename,data_set.data(8737:10416,:));
    elseif i==8
        csvwrite(filename,data_set.data(10225:11904,:));
    elseif i==9
        csvwrite(filename,data_set.data(11713:13344,:));
    elseif i==10
        csvwrite(filename,data_set.data(13153:14832,:));
    elseif i==11
        csvwrite(filename,data_set.data(14641:16272,:));
    elseif i==12
        csvwrite(filename,data_set.data(16081:17760,:));
    end
    end
else
    for i=1:12
        filename=strcat(site,'_',year,'_',num2str(i),'.csv');
    if i==1
        csvwrite(filename,data_set.data(1:1680,:));
    elseif i==2
        csvwrite(filename,data_set.data(1489:3024,:));
    elseif i==3 
        csvwrite(filename,data_set.data(2833:4512,:));
    elseif i==4
        csvwrite(filename,data_set.data(4321:5952,:));
    elseif i==5
        csvwrite(filename,data_set.data(5761:7440,:));
    elseif i==6
        csvwrite(filename,data_set.data(7249:8880,:));
    elseif i==7
        csvwrite(filename,data_set.data(8689:10368,:));
    elseif i==8
        csvwrite(filename,data_set.data(10177:11856,:));
    elseif i==9
        csvwrite(filename,data_set.data(11665:13296,:));
    elseif i==10
        csvwrite(filename,data_set.data(13105:14784,:));
    elseif i==11
        csvwrite(filename,data_set.data(14593:16224,:));
    elseif i==12
        csvwrite(filename,data_set.data(16033:17712,:));
    end
    end
end

    
clear i 
