
mod_data = zeros(16387664,19);
for i = 1104:16388767
    for j = 1:19
        mod_data(i-1103,j) = raw_data(i,j);
    end
end
clear i j 

for i = 1:1103
    for j = 1:19
        mod_data(16387664+i,j) = raw_data(i,j);
    end
end
clear i j 

data = zeros(16388751,19);
for i = 1:16388751
    for j = 1:19
        data(i,j) = mod_data(i,j);
    end
end
clear i j 

temp = 0
for i = 1:16388767
    if mod_data(i,1) > 2013
        temp= temp + 1
    end
end