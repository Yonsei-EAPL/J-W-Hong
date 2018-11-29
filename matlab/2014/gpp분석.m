%%
yr_gb = zeros(1,2); % global annual total (Pg yr-1)
mo_gb = zeros(12,2); % global monthly total (Pg yr-1)
la_gb = zeros(29,2); % latitude total (Pg yr-1)
te_gb = zeros(8,2); % 대륙별 total (Pg yr-1)

%%

gpp = nc_varget('control.GPP.nc','gpp_gb');
gpp_temp = zeros(1,7246);
gpp_m_temp = zeros(12,7246);
gpp(~isfinite(gpp))=0;
for i = 1:2880
    for j =1:7246
        gpp_temp(1,j) = gpp_temp(1,j)+gpp(i,j);
    end
end
clear i j
for i = 1:12
    for j = 1:240
        for k = 1:7246
            gpp_m_temp(i,k) = gpp_m_temp(i,k)+ gpp((i-1)*240+j,k);
        end
    end
end
clear i j k
gpp9999 = zeros(145,192);
gpp9999m = zeros(12,145,192);
for i = 1:7246
    gpp9999(add(1,i),add(2,i)) = gpp_temp(1,i);
    for j =1:12
        gpp9999m(j,add(1,i),add(2,i)) = gpp_m_temp(j,i);
    end
end
clear i j
for i = 1:145
    for j = 1:192
        gpp9999(i,j) = gpp9999(i,j)*area(i,1)*10800*1000;

        for k = 1:12
            gpp9999m(k,i,j) = gpp9999m(k,i,j)*area(i,1)*10800*1000;

        end
    end
end
clear i j k
yr_gb(1,1) = sum(sum(gpp9999))/10^15; % 수정
temp = zeros(145,192);
for i = 1:12
    for j = 1:145
        for k = 1:192
            mo_gb(i,1) = mo_gb(i,1)+gpp9999m(i,j,k); % 수정
            if i==1
                if land(j,k)>0
                    te_gb(land(j,k),1) = te_gb(land(j,k),1) + gpp9999(j,k); % 수정
                end
            end
        end
    end
    mo_gb(i,1) = mo_gb(i,1)/10^15; % 수정
end
clear i j k
for i = 1:29
    for j = 1:5
        for k = 1:192
            la_gb(i,1) = la_gb(i,1) + gpp9999((i-1)*5+j,k); % 수정
        end
    end
end
clear i j k
for i = 1:29
    if i<9
        te_gb(i,1) = te_gb(i,1)/10^15; % 수정
    end
    la_gb(i,1) = la_gb(i,1)/10^15; % 수정
end
clear i
clear gpp gpp9999m gpp_m_temp gpp_temp temp

%% 

resp_p = nc_varget('control.GPP.nc','resp_p_gb');
resp_s = nc_varget('control.GPP.nc','resp_s_gb');
resp = resp_p + resp_s;
clear resp_s resp_p
resp_temp = zeros(1,7246);
resp_m_temp = zeros(12,7246);
resp(~isfinite(resp))=0;
for i = 1:2880
    for j =1:7246
        resp_temp(1,j) = resp_temp(1,j)+resp(i,j);
    end
end
clear i j
for i = 1:12
    for j = 1:240
        for k = 1:7246
            resp_m_temp(i,k) = resp_m_temp(i,k)+ resp((i-1)*240+j,k);
        end
    end
end
clear i j k
resp9999 = zeros(145,192);
resp9999m = zeros(12,145,192);
for i = 1:7246
    resp9999(add(1,i),add(2,i)) = resp_temp(1,i);
    for j =1:12
        resp9999m(j,add(1,i),add(2,i)) = resp_m_temp(j,i);
    end
end
clear i j
for i = 1:145
    for j = 1:192
        resp9999(i,j) = resp9999(i,j)*area(i,1)*10800*1000;
        if resp9999(i,j)>10^15
            resp9999(i,j)=0;
        end
        for k = 1:12
            resp9999m(k,i,j) = resp9999m(k,i,j)*area(i,1)*10800*1000;
            if resp9999m(k,i,j)>10^15
                resp9999m(k,i,j)=0;
            end
        end
    end
end
clear i j k
yr_gb(1,2) = sum(sum(resp9999))/10^15; % 수정
temp = zeros(145,192);
for i = 1:12
    for j = 1:145
        for k = 1:192
            mo_gb(i,2) = mo_gb(i,2)+resp9999m(i,j,k); % 수정
            if i==1
                if land(j,k)>0
                    te_gb(land(j,k),2) = te_gb(land(j,k),2) + resp9999(j,k); % 수정
                end
            end
        end
    end
    mo_gb(i,2) = mo_gb(i,2)/10^15; % 수정
end
clear i j k
for i = 1:29
    for j = 1:5
        for k = 1:192
            la_gb(i,2) = la_gb(i,2) + resp9999((i-1)*5+j,k); % 수정
        end
    end
end
clear i j k
for i = 1:29
    if i<9
        te_gb(i,2) = te_gb(i,2)/10^15; % 수정
    end
    la_gb(i,2) = la_gb(i,2)/10^15; % 수정
end
clear i
clear resp temp resp9999m resp_m_temp resp_temp 
