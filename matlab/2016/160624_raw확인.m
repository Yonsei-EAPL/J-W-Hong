[length a] = size(data);
clear a
temp = 0;
temp2=0;
for i = 1:length
    if i ==1
        temp = temp+1;
        temp2 = data(i,3);
    elseif temp2 ~= data(i,3)
        temp = temp+1;
        temp2 = data(i,3);
    end
        
end
clear i

check = zeros(temp,3);
temp = 0;
temp2=0;
for i = 1:length
    if i ==1
        temp = temp+1;
        temp2 = data(i,3);
        check(temp,1) = temp2;
        check(temp,2) = check(temp,2)+1;
        check(temp,3) = check(temp,3)+1;
    elseif temp2 ~= data(i,3)
        temp = temp+1;
        temp2 = data(i,3);
        check(temp,1) = temp2;
        check(temp,2) = check(temp,2)+1;        
        check(temp,3) = check(temp-1,3)+1;
    else
        check(temp,2) = check(temp,2)+1;
        check(temp,3) = check(temp,3)+1;
    end
        
end
clear i
clear temp temp2