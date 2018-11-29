gpp_raw = gpp;
gpp_raw(~isfinite(gpp_raw))=0;

gpp_post = zeros(1,9323);

for i = 1:28800
    for j = 1:9323
        gpp_post(1,j) = gpp_post(1,j) + gpp_raw(i,j);
    end
end
clear i j

for i = 1:9323
    gpp_post(1,i) = gpp_post(1,i)*3*3600*1000/10;
end
clear i

for i = 1:9323
    gpp_post(1,i) = gpp_post(1,i)*area(loc(1,i),1);
end
clear i

gpp_sum = sum(gpp_post);