n_so2 =0;
n_co =0;
n_no2 =0;
n_pm10 =0;
n_o3 =0;

for i = 1:n_data
    if data2(i,3)>0.15 % so2
        n_so2 = n_so2+1;
    end
    if data2(i,4)>100 % pm10
        n_pm10 = n_pm10+1;
    end
    if data2(i,5)>0.1 % o3
        n_o3 = n_o3+1;
    end
    if data2(i,6)>0.1 % no2
        n_no2 = n_no2+1;
    end
    if data2(i,7)>25 % co
        n_co = n_co+1;
    end    
end
clear i