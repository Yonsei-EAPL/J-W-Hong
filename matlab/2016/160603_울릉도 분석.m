co = raw;
co = sortrows(co,3);
co = sortrows(co,4);
co(1:7624,:) = [];
% length = 113952;
% for i = 1:length
%     if (co(i,3)==0)&&(co(i,4)==0)
%         co(i,:)=[];
%     end
%     i=i-1;
% end
% clear i

co2 =zeros(106328,3);
for i = 1:106328
    co2(i,1) = co(i,2);%(co(i,1)-2001)*12+co(i,2);
    co2(i,2) = co(i,3);
    co2(i,3) = co(i,4);
end
clear i
co2 = sortrows(co2,1);
co2(106328,:)=[];

n_data = zeros(156,4);
for i = 1:156
    n_data(i,1) = i;
end
clear i
for i = 1:106327
    n_data(co2(i,1),2) = n_data(co2(i,1),2)+1;
end 
clear i

for i = 1:12%156
    if i==1
        [a0 a1] = linear_regression(co2(1:n_data(i,2),2),co2(1:n_data(i,2),3));
        a = corrcoef(co2(1:n_data(i,2),2),co2(1:n_data(i,2),3));
        n_data(i,3) = a1;
        n_data(i,4) = a(1,2)^2;
    else
        [a0 a1] = linear_regression(co2(sum(n_data(1:i-1,2)):sum(n_data(1:i,2)),2),co2(sum(n_data(1:i-1,2)):sum(n_data(1:i,2)),3));
        a = corrcoef(co2(sum(n_data(1:i-1,2)):sum(n_data(1:i,2)),2),co2(sum(n_data(1:i-1,2)):sum(n_data(1:i,2)),3));
        n_data(i,3) = a1;
        n_data(i,4) = a(1,2)^2;
    end
end
clear i

