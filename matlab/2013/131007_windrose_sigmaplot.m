[size_n size_v] = size(data);
result = zeros(3600,2);

wd_max = data(size_n,2);
wd_min = data(1,1);

if wd_min - (360 - wd_max) > 0
    temp0 = data(size_n,3);
    for i = 1:3600
       x = (i-1)*0.1;
       result(i,1) = x;
       if i ==1
           y = temp0;
           temp = 0;
       else
           if temp < size_n
               if data(temp+1,1)<x
                   y = data(temp+1,3);
                   temp = temp+1;
               end
           else
               if data(temp,2)<x
                   y = temp0;
               end
           end
           
       end
       result(i,2) = y;
    end
else
    temp0 = data(1,3);
    for i = 1:3600
       x = (i-1)*0.1;
       result(i,1) = x;
       if i ==1
           y = temp0;
           temp = 0;
       else
           if temp < size_n
               if data(temp+1,1)<x
                   y = data(temp+1,3);
                   temp = temp+1;
               end
           else
               if data(temp,2)<x
                   y = temp0;
               end
           end
           
       end
       result(i,2) = y;
    end    
end
clear x y i temp temp0 wd_max wd_min size_n size_v;
