fp_coab = fp_co;

for i = 1:21640
    for j = 1:43
        fp_coab(i,j) = abs(fp_coab(i,j));
    end
end
clear i j

for i = 1:21640
    temp = sum(fp_coab(i,:));
    for j = 1:43
        fp_coab(i,j) = fp_coab(i,j)/temp;
    end
end
clear i j