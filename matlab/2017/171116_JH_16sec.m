
data = importdata('BS_lvl4_2016.txt');
data=strrep(data,'-','');
data=strrep(data,':','');
data=strrep(data,'  ',',');
data=strrep(data,',,',',');
data=strrep(data,',,',',');
data=strrep(data,',,',',');
data=strrep(data,',,',',');
data=strrep(data,',,',',');
data=strrep(data,',,',',');
data=strrep(data,',,',',');
data=strrep(data,',,',',');
data=strrep(data,' ','');

data=regexp(data,',','split');
result=zeros(length(data),3);
for i=1:length(data)
    if size(data{i,:},2)~=5
    else
        temp=data{i,1}(1,1);
        result(i,1)=str2num(temp{1,1});
        temp=data{i,1}(1,2);
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
end
clear i temp2 ss hh mm

result2 = result(:,3);
result2 = unique(result2);
result2(:,2) = 0;
result2(:,3) = 0;
for i = 1:length(result)
    if (result(i,2)~=9999999)&&(result(i,2)~=0)
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



