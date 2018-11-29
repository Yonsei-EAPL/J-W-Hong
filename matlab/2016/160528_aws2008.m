temp = 0;

for i = 1:5596499
    if aws2008(i,2)==371
        temp = temp+1;
    end
end
clear i

temp=0;
d371 = zeros(3330,11);
for i = 1:5596499
    if aws2008(i,2)==371
       temp = temp+1;
       for j = 1:11
        d371(temp,j) = aws2008(i,j);
       end
    end
end
clear i