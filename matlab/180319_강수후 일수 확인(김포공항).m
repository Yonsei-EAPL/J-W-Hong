x(:,3) = 0;

for i = 1:length(x)
    if x(i,2)>0
        x(i,3)=0;
    else
        x(i,3) = x(i-1,3)+1;
    end
end
clear i 