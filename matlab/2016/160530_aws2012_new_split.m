%% with 2012_new

stn_n = 418;

temp=0;
for i = 1:3977018
    if aws2012(i,2)==stn_n
        temp = temp+1;
    end
end
clear i

data = zeros(temp,11);
temp=0;
for i = 1:3977018
    if aws2012(i,2)==stn_n
        temp = temp+1;
        for j = 1:11
            data(temp,j) = aws2012(i,j+1);
        end
    end
    clear j
end
clear i temp

