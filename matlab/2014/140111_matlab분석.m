YR_iso = zeros(360,720);
JJA_iso = zeros(360,720);
DJF_iso = zeros(360,720);

for i = 1:360
    for j = 1:720
        YR_iso(i,j) = isoprene2004(13,i,j)/10^6*area(i,1)*3;
        for k= 1:3
            JJA_iso(i,j) = JJA_iso(i,j) + isoprene2004(5+k,i,j)/10^6*area(i,1)*3;
            if k==1
                DJF_iso(i,j) = DJF_iso(i,j) + isoprene2004(12,i,j)/10^6*area(i,1)*3;
            else
                DJF_iso(i,j) = DJF_iso(i,j) + isoprene2004(k-1,i,j)/10^6*area(i,1)*3;
            end
        end
    end
end
clear i j k

a = zeros(360,720);
monthly_iso = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = isoprene2004(k,i,j)*area(i,1)/10^6*3;
        end
    end
    clear i j
    monthly_iso(k,1) = sum(sum(a));
end
clear a k