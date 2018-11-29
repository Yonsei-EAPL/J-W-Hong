min_T=min(a(:,1));
max_T=max(a(:,1));
ra_T = (max_T - min_T)/9;

min_W=min(a(:,2));
max_W=max(a(:,2)); 
ra_W = (max_W - min_W)/9;

result = zeros(10,10);
result_n = zeros(10,10);
[x y] = size(a);

for i = 1:x
    temp_x = (a(i,1) - min_T)/ra_T;
    temp_x = temp_x - mod(temp_x,1);
    if temp_x ==0
        temp_x = 1;
    elseif temp_x >10
        temp_x = 10;
    end
    
    temp_y = (a(i,2) - min_W)/ra_W;
    temp_y = temp_y - mod(temp_y,1);
    if temp_y ==0
        temp_y = 1;
    elseif temp_y >10
        temp_y = 10;
    end
    
    result(temp_x,temp_y) =result(temp_x,temp_y) + a(i,3);
    result_n(temp_x,temp_y) =result_n(temp_x,temp_y) + 1;
end
clear i temp_x temp_y ra_T ra_W x y

for i = 1:10
    for j = 1:10
        if result_n(i,j)>0
            result(i,j) = result(i,j)/result_n(i,j);
        end
    end
end
clear i j
clear min_W max_W min_T max_T