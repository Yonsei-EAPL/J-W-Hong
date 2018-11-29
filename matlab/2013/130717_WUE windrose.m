result = zeros(36000,2);
for i = 1:36000
    result(i,1) = i/100;
end

temp = zeros(44,1);
for i = 1:22
    for j = 1:2
        temp((i-1)*2+j,1)= data(i,j);
    end
end

temp2 = zeros(23,2);
for i = 2:22
    temp2(i,1) = (data(i-1,2)+data(i,1))/2;
    temp2(i,2) = data(i-1,3);
end
temp2(1,1) = 0;
temp2(1,2) = data(1,3);
temp2(23,1)= 360;
temp2(23,2)= data(22,3);

for i = 1:36000
    temp=0;
    for j = 1:23
        if (temp2(j,1)<(i/100))&&(temp2(j+1,1)>(i/100))
            temp = temp2(j+1,2);
        end
    end
    result(i,2) = temp;
end

clear i j 

