n = 10511; % length

% x : raw data

y = zeros(n,2); % 0.005 interval
for i = 1:n
    y(i,1) = x(i,1)-mod(x(i,1),0.005)+0.005;
    y(i,2) = x(i,2)-mod(x(i,2),0.005)+0.005;
end
clear i

xx= zeros(200,5);
for i = 1:200
    xx(i,1) = i*0.005;
end
clear i

for i = 1:n
    temp = int64(y(i,1)/0.005);
    xx(temp,2) = xx(temp,2)+1;
    temp = int64(y(i,2)/0.005);
    xx(temp,3) = xx(temp,3)+1;
end
clear i

for i = 1:200
    if i ==1
        xx(i,4) = xx(i,2)/10511;
        xx(i,5) = xx(i,3)/10511;
    else
        xx(i,4) = xx(i-1,4) + xx(i,2)/10511;
        xx(i,5) = xx(i-1,5) + xx(i,3)/10511;
    end
end
clear i

xxx= zeros(100,100);
for i = 1:n
    temp1 = int64(y(i,1)/0.005);
    temp2 = int64(y(i,2)/0.005);
    if mod(temp1,2)==1
        temp1 = (temp1+1)/2;
    else
        temp1 = temp1/2;
    end
    if mod(temp2,2)==1
        temp2 = (temp2+1)/2;
    else
        temp2 = temp2/2;
    end    
    xxx(temp1,temp2) = xxx(temp1,temp2) +1;
end
clear i
xxx = xxx/n*100;

