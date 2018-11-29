folderDir = 'H:\b1_EAPL\b1_JW_Observation\b1_obs_보성\보성운고계(CL51)\'; 
folderName = dir(fullfile(folderDir))';
folderlen = length(folderName);
result_all = 0;
for tt=3:4%folderlen
    dataName = dir(fullfile(folderDir, folderName(tt).name,'CEILOMETER_1_LEVEL_3_DEFAULT_*'))'; 
    total_result=cell(length(dataName),2);    
    for t=1:length(dataName)
        data = importdata(fullfile(folderDir, folderName(tt).name, dataName(t).name));
        data(1,:)=[];
        data(1,:)=[];

        data=strrep(data,'-','');
        data=strrep(data,':','');
        data=strrep(data,' ','');
        data=strrep(data,'DEVICE_1,','');
        data=strrep(data,'_','');
        data=regexp(data,',','split');
        result=zeros(length(data),3);

        for i=1:length(data)
            temp=data{i,1}(1,1);
            result(i,1)=str2num(temp{1,1});
            temp=data{i,1}(1,5);
            result(i,2)=str2num(temp{1,1});
            clear temp

            ss = mod(result(i,1),100);
            mm = mod(result(i,1),10^4);
            mm = (mm - ss)/100;
            hh = mod(result(i,1),10^6);
            hh = (hh - mm*10^4 - ss*100)/10^4;

            temp2=0;
            if (ss==0)&&((mm==0)||(mm==30))
                if mm==0
                    if hh==0
                        mm=15;
                    elseif hh==23
                        mm=45;
                    else
                        hh=hh-1;
                        mm=45;
                        temp2=1;
                    end

                else
                    mm=15;
                end
            else
                if mm==0
                    mm=15;
                elseif mm==30
                    mm=45;
                elseif mm<30
                    mm=15;
                else
                    mm=45;
                end
            end
            if temp2==1
                result(i,3) = result(i,1)-mod(result(i,1),10^6);
                result(i,3) = result(i,3)+hh*10^4+mm*100;            
            else
                result(i,3) = result(i,1)-mod(result(i,1),10^4);
                result(i,3) = result(i,3)+mm*100;
            end
        end
        clear i temp2 ss hh mm

        result2 = result(:,3);
        result2 = unique(result2);
        result2(:,2) = 0;
        result2(:,3) = 0;

        for i = 1:length(result)
            if result(i,2)~=999
                xx = find(result2==result(i,3));
                result2(xx,2) = result2(xx,2)+result(i,2);
                result2(xx,3) = result2(xx,3)+1;
            end
        end
        clear i xx
        for i = 1:length(result2)
            if result2(i,3)==0
                result2(i,2)=-999;
            else
                result2(i,2)=result2(i,2)/result2(i,3);
            end
        end
        clear i result
        total_result{t,1}= dataName(t).name;
        total_result{t,2}= result2;
        clear result2
        clear data
    end
    clear t

    result = 0;
    for i = 1:length(total_result)
        if i==1
            result = total_result{i,2};
        else
            result = [result; total_result{i,2}];
        end
    end
    clear i 
    clear total_result dataDir dataName
    result_all = [result_all; result];
    clear result
end
clear tt folderDir folderlen folderName



