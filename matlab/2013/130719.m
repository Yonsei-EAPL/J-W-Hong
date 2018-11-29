%% for windrose_sigmaplot

[size_x size_y] = size(data);

result = zeros(3600,2);
for i = 1:3600
    result(i,1) = i/10;
end

temp = zeros(size_x+2,2);
temp(1,2) = data(1,3);

last = 360 + (data(1,1)+data(size_x,2)-360)/2;
temp(size_x+1,1) = last;
clear last
temp(size_x+1,2) = data(size_x,3);

temp(size_x+2,1) = 360;
temp(size_x+2,2) = data(1,3);

for i = 1:size_x-1
    temp(i+1,1) = (data(i,2)+data(i+1,1))/2;
    temp(i+1,2) = data(i,3);
end

for i = 1:3600
    temp2=0;
    for j = 1:size_x
        if (i/10)>temp(j,1)&&(i/10)<temp(j+1,1)
            temp2 = temp(j+1,2);
        end
    end
    if (i/10)>temp(size_x+1,1)
        temp2 = temp(size_x+2,2);
    end
    result(i,2) = temp2;
end
clear i j    
    

%%

for i = 1:3600
    result(i,2) = abs(result(i,2));
end
clear i