[n, x] = hist(data,100);
for i = 1:100
    n(1,i) = n(1,i)*100/43769;
end
clear i
bar(x,n)