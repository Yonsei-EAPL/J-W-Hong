[size_n size_var] = size(data);

temp = 0;
for i= 1:size_n
    if (data(i,1)==0)&&(temp==0)
        temp = 1;
    elseif (data(i,1)==0)&&(temp>0)
        temp = temp +1;
    elseif (data(i,1)>0)&&(temp==0)
        temp2 = data(i,1);
    elseif (data(i,1)>0)&&(temp>0)
        temp3 = data(i,1);
        for j = 1:temp
            slope = (temp3-temp2)/(temp+1);
            data(i-(temp+1)+j,1) = temp2 + slope*j;
        end
        temp = 0;
        temp2 = temp3;
        temp3 = 0;
    end
end

clear i j temp temp2 temp3