% xx = zeros(163632,3);
% 
% for i = 1:163631
%     xx(i+1,1) = x(i,1);
%     xx(i+1,2) = x(i,2);
%     xx(i+1,3) = x(i,3);
% end
% clear i
% 
% x=xx;
% clear xx

result = zeros(20454,3);

for i = 1:20454
    temp = zeros(8,1);
    for j=1:8
        temp(j,1) = x((i-1)*8+j,3);
    end
    clear j
    [tt yy] = max(temp);
    
    result(i,1) = x((i-1)*8+1,1);
    result(i,2) = tt;
    result(i,3) = (yy-1)*3;
    clear tt yy
    
end
clear i 