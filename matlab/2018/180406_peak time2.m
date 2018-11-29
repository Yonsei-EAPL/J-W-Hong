xx= zeros(20455,4);

temp = 0;
temp2 = 0;

for i = 1:length(x(:,1))
    if i==1
        temp = temp+1;
        temp2 = x(1,1);
        xx(temp,1) = x(i,1);
        xx(temp,2) = xx(temp,2)+1;
        xx(temp,3) = x(i,2);
        xx(temp,4) = x(i,3);
    else
        if temp2 ~=x(i,1)
            temp = temp+1;
            temp2 = x(i,1);
            xx(temp,1) = x(i,1);
            xx(temp,2) = xx(temp,2)+1;
            xx(temp,3) = x(i,2);
            xx(temp,4) = x(i,3);
        else
            xx(temp,2) = xx(temp,2)+1;            
            if xx(temp,4)<x(i,3)
                xx(temp,3) = x(i,2);
                xx(temp,4) = x(i,3);
            end
        end
    end
end
clear i