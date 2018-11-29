%% for 100

result = zeros(10000,2);
temp =0;
for i = 1:10000
    temp = 0;
    for j = 1:1000
        if (y4(i,j)<100)&&(temp==0)
            result(i,1) = log10(j*100);
            result(i,2) = log10(i*10000);
            temp=temp+1;
        end        
    end
end
clear i j
    