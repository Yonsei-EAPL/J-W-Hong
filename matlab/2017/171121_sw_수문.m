nn = 236;
data = ascat;

%nn=10324;
%data = crs(:,1);
%data = crs(:,2);
% data(:,2)=0;

x= zeros(20,4);
for i = 1:20
    x(i,1) = (i-1)*0.05+0.025;
end
clear i

for i = 1:nn
   if data(i,1)==1
       data(i,2) = 20;
   else
       data(i,2) = ceil((data(i,1) - mod(data(i,1),0.05)+0.05)/0.05);
   end
end
clear i

for i = 1:nn
   x(data(i,2),2) =  x(data(i,2),2)+1;
end
clear i

for i = 1:20
    x(i,3) = x(i,2)/nn *100;
    if i==1
        x(i,4) = x(i,3)/100;
    else
        x(i,4) = x(i-1,4)+x(i,3)/100;
    end
end
clear i 
