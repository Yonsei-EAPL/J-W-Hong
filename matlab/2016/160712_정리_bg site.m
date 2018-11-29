% 10007 - pm10 ug m-3
% 10006 - o3 ppm 0.039
% 10003 - NO2 ppm 0.034
% 10002 - CO ppm 0.5
% 10001 - SO ppm 0.005
 
data2=sortrows(data,[2 1]);

[a b] = size(data2);
clear b

temp = 0;
temp2=0;
for i = 1:a
    if i==1
        temp = temp+1;
        temp2 = data2(i,2);
    else
        if temp2 ~=data2(i,2)
            temp=temp+1;
            temp2=data2(i,2);
        end
    end
end
clear i

data_inf = zeros(temp,6);
temp = 0;
temp2=0;
for i = 1:a
    if i==1
        temp = temp+1;
        data_inf(temp,1) = data2(i,2);
        temp2 = data2(i,2);
    else
        if temp2 ~=data2(i,2)
            temp=temp+1;
            data_inf(temp,1) = data2(i,2);
            temp2=data2(i,2);
        end
    end
end
clear i

temp = 0;
temp2=0;
for i = 1:a
    if i==1
        temp = temp+1;
        temp2 = data2(i,2);
    else
        if temp2 ~=data2(i,2)
            temp=temp+1;
            temp2=data2(i,2);
        end
    end
    if data2(i,3)==10001
        data_inf(temp,2) = data_inf(temp,2)+1;
    elseif data2(i,3)==10002
        data_inf(temp,3) = data_inf(temp,3)+1;
    elseif data2(i,3)==10003
        data_inf(temp,4) = data_inf(temp,4)+1;
    elseif data2(i,3)==10006
        data_inf(temp,5) = data_inf(temp,5)+1;
    else
        data_inf(temp,6) = data_inf(temp,6)+1;        
    end
end
clear i
clear temp temp2

ans = sum(data_inf);
aa=ans(1,2);
clear ans

data3 = zeros(aa,11);

for i=1:aa*5
    if i==1
        temp=1;
        temp2=data2(i,1);
        data3(1,1) = data2(1,1);
        data3(1,2) = data2(1,2);
    else
        if temp2~=data2(i,1)
            temp=temp+1;
            temp2=data2(i,1);
            data3(temp,1) = data2(i,1);
            data3(temp,2) = data2(i,2);
        end
    end
    if data2(i,3)==10001
        data3(temp,3)=data2(i,4);
    elseif data2(i,3)==10002
        data3(temp,4)=data2(i,4);
    elseif data2(i,3)==10003
        data3(temp,5)=data2(i,4);
    elseif data2(i,3)==10006
        data3(temp,6)=data2(i,4);
    else
        data3(temp,7)=data2(i,4);
    end   
    data3(temp,11) = mod(data3(temp,1),100);
    data3(temp,10) = mod(data3(temp,1),10000);
    data3(temp,10) = (data3(temp,10)-mod(data3(temp,10),100))/100;
    data3(temp,9) = mod(data3(temp,1),1000000);
    data3(temp,9) = (data3(temp,9)-mod(data3(temp,9),10000))/10000; 
    data3(temp,8) = mod(data3(temp,1),100000000);
    data3(temp,8) = (data3(temp,8)-mod(data3(temp,8),1000000))/1000000;     
   
end
clear i
clear temp temp2

for i=1:aa
    data3(i,2) = data3(i,2)-831400;
    data3(i,1) = (data3(i,1)-mod(data3(i,1),100000000))/100000000;
end
clear i


csvwrite('data_bg.dat', data3);



% clear a