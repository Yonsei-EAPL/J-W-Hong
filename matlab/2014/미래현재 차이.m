iso_rcp = zeros(145,192);
iso_hist = zeros(145,192);
ter_rcp = zeros(145,192);
ter_hist = zeros(145,192);
iso_diff = zeros(145,192);
ter_diff = zeros(145,192);

for i = 1:145
    for j = 1:192
        iso_rcp(i,j)=(iso2090(i,j)+iso2091(i,j)+iso2092(i,j)+iso2093(i,j)+iso2094(i,j)+iso2095(i,j)+iso2096(i,j)+iso2097(i,j)+iso2098(i,j)+iso2099(i,j))/10;
        iso_hist(i,j)=(iso1996(i,j)+iso1997(i,j)+iso1998(i,j)+iso1999(i,j)+iso2000(i,j)+iso2001(i,j)+iso2002(i,j)+iso2003(i,j)+iso2004(i,j)+iso2005(i,j))/10; 
        iso_diff(i,j) = iso_rcp(i,j)-iso_hist(i,j);
        ter_rcp(i,j)=(ter2090(i,j)+ter2091(i,j)+ter2092(i,j)+ter2093(i,j)+ter2094(i,j)+ter2095(i,j)+ter2096(i,j)+ter2097(i,j)+ter2098(i,j)+ter2099(i,j))/10;
        ter_hist(i,j)=(ter1996(i,j)+ter1997(i,j)+ter1998(i,j)+ter1999(i,j)+ter2000(i,j)+ter2001(i,j)+ter2002(i,j)+ter2003(i,j)+ter2004(i,j)+ter2005(i,j))/10;         
        ter_diff(i,j) = ter_rcp(i,j)-ter_hist(i,j);
    end
end
clear i j

for i = 1:145
    for j = 1:192
        if land(i,j)==0
            iso_diff(i,j) = -10^(12);
            ter_diff(i,j) = -10^(11);
        end
    end
end
clear i j