for i = 1:4175
    j=find(x(:,1)==y(i,1));
    x(j,2) = y(i,2);
end
clear i

for i = 1:20638
    if x(i,2)==0
        x(i,2)=-999;
    end
end 
clear i