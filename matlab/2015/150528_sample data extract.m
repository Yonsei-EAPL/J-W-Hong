a=raw1;%(:,10);
b=raw2;%(:,10);
clear raw1 raw2
co2 = zeros(length(a)+length(b),3);

for i = 1:length(a)
    co2(i,1) = a(i,1);
end
for i = 1:length(b)
    co2(i+length(a),1) = b(i,1);
end
clear a b
clear i