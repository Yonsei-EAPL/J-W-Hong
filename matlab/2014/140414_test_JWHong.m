temp = zeros(16008001,1);

for i = 1:4001
    for j = 1:4001
        temp((i-1)*4001+j,1) = east(i,j);
    end
end
clear i j

y = quantile(temp,[.025 .25 .50 .75 .975]);