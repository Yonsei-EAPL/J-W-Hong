%% using control_result (EAP class material, sc10, 2011)

a = zeros(360,720);
monthly_gpp = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_gpp(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_gpp(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_iso = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_iso(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_iso(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_re = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_re(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_re(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_re_p = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_re_p(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_re_p(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_re_s = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_re_s(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_re_s(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_npp = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_npp(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_npp(k,1) = sum(sum(a));
end
clear a k

YR_gpp = zeros(360,720);
JJA_gpp = zeros(360,720);
DJF_gpp = zeros(360,720);
for i = 1:360
    for j = 1:720
        YR_gpp(i,j) = con_gpp(13,i,j);
        for k = 1:3
            JJA_gpp(i,j) = JJA_gpp(i,j) +con_gpp(5+k,i,j);
            if k==1
                DJF_gpp(i,j) = DJF_gpp(i,j) +con_gpp(12,i,j);
            else
                DJF_gpp(i,j) = DJF_gpp(i,j) +con_gpp(k-1,i,j);
            end
        end
    end
end
clear i j k

YR_re_s = zeros(360,720);
JJA_re_s = zeros(360,720);
DJF_re_s = zeros(360,720);
for i = 1:360
    for j = 1:720
        YR_re_s(i,j) = con_re_s(13,i,j);
        for k = 1:3
            JJA_re_s(i,j) = JJA_re_s(i,j) +con_re_s(5+k,i,j);
            if k==1
                DJF_re_s(i,j) = DJF_re_s(i,j) +con_re_s(12,i,j);
            else
                DJF_re_s(i,j) = DJF_re_s(i,j) +con_re_s(k-1,i,j);
            end
        end
    end
end
clear i j k

YR_iso = zeros(360,720);
JJA_iso = zeros(360,720);
DJF_iso = zeros(360,720);
for i = 1:360
    for j = 1:720
        YR_iso(i,j) = con_iso(13,i,j);
        for k = 1:3
            JJA_iso(i,j) = JJA_iso(i,j) +con_iso(5+k,i,j);
            if k==1
                DJF_iso(i,j) = DJF_iso(i,j) +con_iso(12,i,j);
            else
                DJF_iso(i,j) = DJF_iso(i,j) +con_iso(k-1,i,j);
            end
        end
    end
end
clear i j k


contourf(YR_gpp,100,'LineStyle','none')
hold on
plot(m_lon(:,1),m_lat(:,1),'k','LineWidth',0.5)

contourf(YR_gpp,100,'LineStyle','none')
contourf(JJA_gpp,100,'LineStyle','none')
contourf(DJF_gpp,100,'LineStyle','none')

contourf(YR_iso,100,'LineStyle','none')
contourf(JJA_iso,100,'LineStyle','none')
contourf(DJF_iso,100,'LineStyle','none')

contourf(YR_re_s,100,'LineStyle','none')
contourf(JJA_re_s,100,'LineStyle','none')
contourf(DJF_re_s,100,'LineStyle','none')


%% using control (sc 5, BVOC project, 2011)
% Import the Result
gpp = nc_varget('control.monthly.nc','gpp_gb');
npp = nc_varget('control.monthly.nc','npp_gb');
respp = nc_varget('control.monthly.nc','resp_p_gb');
resps = nc_varget('control.monthly.nc','resp_s_gb');
isoprene = nc_varget('control.monthly.nc','isoprene');

% Convert to Monthly Accumulated Result
for i = 1:12
    for j = 1:67209
        gpp(i,j) = gpp(i,j)*DOM(i,1)*(24*3600*10^3); % kgCm-2s-1 ; gCm-2month-1
        npp(i,j) = npp(i,j)*DOM(i,1)*(24*3600*10^3); % kgCm-2s-1 ; gCm-2month-1
        respp(i,j) = respp(i,j)*DOM(i,1)*(24*3600*10^3); % kgCm-2s-1 ; gCm-2month-1
        resps(i,j) = resps(i,j)*DOM(i,1)*(24*3600*10^3); % kgCm-2s-1 ; gCm-2month-1
        isoprene(i,j) = isoprene(i,j)*DOM(i,1)*(24*3600*10^3); % kgCm-2s-1 ; gCm-2month-1
    end
end
clear i j

% Make Grid 
gpp_g = zeros(13,360,720); % 13 means annual sum
npp_g = zeros(13,360,720);
respp_g = zeros(13,360,720);
resps_g = zeros(13,360,720);
isoprene_g = zeros(13,360,720);

for i = 1:12
    for j = 1:67209
        gpp_g(i,add_result(j,1),add_result(j,2))=gpp(i,j);
        gpp_g(13,add_result(j,1),add_result(j,2))=gpp_g(13,add_result(j,1),add_result(j,2))+gpp(i,j);
        npp_g(i,add_result(j,1),add_result(j,2))=npp(i,j);
        npp_g(13,add_result(j,1),add_result(j,2))=npp_g(13,add_result(j,1),add_result(j,2))+npp(i,j);
        respp_g(i,add_result(j,1),add_result(j,2))=respp(i,j);
        respp_g(13,add_result(j,1),add_result(j,2))=respp_g(13,add_result(j,1),add_result(j,2))+respp(i,j);
        resps_g(i,add_result(j,1),add_result(j,2))=resps(i,j);
        resps_g(13,add_result(j,1),add_result(j,2))=resps_g(13,add_result(j,1),add_result(j,2))+resps(i,j);
        isoprene_g(i,add_result(j,1),add_result(j,2))=isoprene(i,j);
        isoprene_g(13,add_result(j,1),add_result(j,2))=isoprene_g(13,add_result(j,1),add_result(j,2))+isoprene(i,j);
    end
end
clear i j
clear gpp npp respp resps isoprene
gpp = gpp_g;
npp = npp_g;
respp = respp_g;
resps = resps_g;
isoprene = isoprene_g;

% Match with Map (Turn 360-grid)
for k = 1:13
    for i = 1:360
        for j = 1:360
            gpp(k,i,j) = gpp_g(k,i,j+360);
            npp(k,i,j) = npp_g(k,i,j+360);
            respp(k,i,j) = respp_g(k,i,j+360);
            resps(k,i,j) = resps_g(k,i,j+360);
            isoprene(k,i,j) = isoprene_g(k,i,j+360);
        end
    end
end
for k = 1:13
    for i = 1:360
        for j = 361:720
            gpp(k,i,j) = gpp_g(k,i,j-360);
            npp(k,i,j) = npp_g(k,i,j-360);
            respp(k,i,j) = respp_g(k,i,j-360);
            resps(k,i,j) = resps_g(k,i,j-360);
            isoprene(k,i,j) = isoprene_g(k,i,j-360);
        end
    end
end
clear i j k
clear gpp_g npp_g respp_g resps_g isoprene_g

resp = respp + resps;

% Final Rename
con_gpp = gpp;
clear gpp
con_npp = npp;
clear npp
con_re = resp;
clear resp
con_re_p = respp;
clear respp
con_re_s = resps;
clear resps
con_iso = isoprene;
clear isoprene

a = zeros(360,720);
monthly_gpp = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_gpp(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_gpp(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_iso = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_iso(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_iso(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_re = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_re(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_re(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_re_p = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_re_p(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_re_p(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_re_s = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_re_s(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_re_s(k,1) = sum(sum(a));
end
clear a k

a = zeros(360,720);
monthly_npp = zeros(13,1);
for k = 1:13
    for i = 1:360
        for j = 1:720
            a(i,j) = con_npp(k,i,j)*area(i,1);
        end
    end
    clear i j
    monthly_npp(k,1) = sum(sum(a));
end
clear a k

YR_gpp = zeros(360,720);
JJA_gpp = zeros(360,720);
DJF_gpp = zeros(360,720);
for i = 1:360
    for j = 1:720
        YR_gpp(i,j) = con_gpp(13,i,j);
        for k = 1:3
            JJA_gpp(i,j) = JJA_gpp(i,j) +con_gpp(5+k,i,j);
            if k==1
                DJF_gpp(i,j) = DJF_gpp(i,j) +con_gpp(12,i,j);
            else
                DJF_gpp(i,j) = DJF_gpp(i,j) +con_gpp(k-1,i,j);
            end
        end
    end
end
clear i j k

YR_re_s = zeros(360,720);
JJA_re_s = zeros(360,720);
DJF_re_s = zeros(360,720);
for i = 1:360
    for j = 1:720
        YR_re_s(i,j) = con_re_s(13,i,j);
        for k = 1:3
            JJA_re_s(i,j) = JJA_re_s(i,j) +con_re_s(5+k,i,j);
            if k==1
                DJF_re_s(i,j) = DJF_re_s(i,j) +con_re_s(12,i,j);
            else
                DJF_re_s(i,j) = DJF_re_s(i,j) +con_re_s(k-1,i,j);
            end
        end
    end
end
clear i j k

YR_iso = zeros(360,720);
JJA_iso = zeros(360,720);
DJF_iso = zeros(360,720);
for i = 1:360
    for j = 1:720
        YR_iso(i,j) = con_iso(13,i,j);
        for k = 1:3
            JJA_iso(i,j) = JJA_iso(i,j) +con_iso(5+k,i,j);
            if k==1
                DJF_iso(i,j) = DJF_iso(i,j) +con_iso(12,i,j);
            else
                DJF_iso(i,j) = DJF_iso(i,j) +con_iso(k-1,i,j);
            end
        end
    end
end
clear i j k











