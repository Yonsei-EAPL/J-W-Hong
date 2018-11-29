
x(1,2) = 0;
x(1,3) = 0;
for i = 1:length(x(:,1))
    if x(i,1)<=-12
        x(i,2)=1;
        if x(i-1,2)==1
            x(i,3)=1;
        end
    end
    if i>1
        if (((x(i-1,1)-x(i,1))>=10)&&(x(i,1)<=3))
            x(i,2)=1;
            x(i,3)=1;
        end
    end
end
clear i