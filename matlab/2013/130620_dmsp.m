%% SEOUL separate
display('SEOUL separate')

size_si = 200;
size_sj = 200;

seoul = zeros(31,size_si,size_sj);
data_mask_seoul = zeros(size_si,size_sj);

for i = 1:size_si
    for j = 1:size_sj
        seoul(1,i,j) = f101992(i+775+51,j+540+181);
        seoul(2,i,j) = f101993(i+775+51,j+540+181);        
        seoul(3,i,j) = f101994(i+775+51,j+540+181);
        seoul(4,i,j) = f121994(i+775+51,j+540+181);
        seoul(5,i,j) = f121995(i+775+51,j+540+181);
        seoul(6,i,j) = f121996(i+775+51,j+540+181);        
        seoul(7,i,j) = f121997(i+775+51,j+540+181);
        seoul(8,i,j) = f121998(i+775+51,j+540+181);
        seoul(9,i,j) = f121999(i+775+51,j+540+181);
        seoul(10,i,j) = f141997(i+775+51,j+540+181);
        seoul(11,i,j) = f141998(i+775+51,j+540+181);
        seoul(12,i,j) = f141999(i+775+51,j+540+181);
        seoul(13,i,j) = f142000(i+775+51,j+540+181);
        seoul(14,i,j) = f142001(i+775+51,j+540+181);
        seoul(15,i,j) = f142002(i+775+51,j+540+181);
        seoul(16,i,j) = f142003(i+775+51,j+540+181);
        seoul(17,i,j) = f152000(i+775+51,j+540+181);
        seoul(18,i,j) = f152001(i+775+51,j+540+181);
        seoul(19,i,j) = f152002(i+775+51,j+540+181);
        seoul(20,i,j) = f152003(i+775+51,j+540+181);
        seoul(21,i,j) = f152004(i+775+51,j+540+181);
        seoul(22,i,j) = f152005(i+775+51,j+540+181);
        seoul(23,i,j) = f152006(i+775+51,j+540+181);
        seoul(24,i,j) = f152007(i+775+51,j+540+181);
        seoul(25,i,j) = f162004(i+775+51,j+540+181);
        seoul(26,i,j) = f162005(i+775+51,j+540+181);
        seoul(27,i,j) = f162006(i+775+51,j+540+181);
        seoul(28,i,j) = f162007(i+775+51,j+540+181);
        seoul(29,i,j) = f162008(i+775+51,j+540+181);
        seoul(30,i,j) = f162009(i+775+51,j+540+181);
        seoul(31,i,j) = f182010(i+775+51,j+540+181);
        data_mask_seoul(i,j) = data_mask(i+775+51,j+540+181);
    end
end
clear i j 


%% percentile
display('percentile')

for i = 1:31
    for j = 1:size_si
        for k = 1:size_sj
            seoul(i,j,k) = (seoul(i,j,k)/63)*100;
        end
    end
end
clear i j k


%% threshold for each satellite
figure(1)
for i = 1:3
    ii = i;
    if i == 1
        i = 30;
    elseif i == 2
        i = 1;
    elseif i == 3
        i = 24;
    end
    
    temp = zeros(size_si,size_sj);
    th = zeros(100,1);
    for j = 1:100
        for k = 1:size_si
            for l = 1:size_sj
                temp(k,l) = seoul(i,k,l);
            end
        end
        for k = 1:size_si
            for l = 1:size_sj
                if temp(k,l)<j
                    temp(k,l) =0;
                end
            end
        end
        peri = bwperim(temp);
        th(j,1) = sum(sum(peri));
    end
    a=max(th);
    for j = 1:100
        th(j,1) = th(j,1)/a;
    end
    
    plot(th(:,1))
    hold on
    
    i = ii;
end
clear i j k a temp th peri ii l


%% pictures

%integrate data set
total = zeros(31,size_i, size_j);
for i = 1:31
    for j = 1:size_i
        for k = 1:size_j
            if i == 1
                total(i,j,k) = f101992_land(j,k);
            elseif i ==2
                total(i,j,k) = f101993_land(j,k);            
            elseif i ==3
                total(i,j,k) = f101994_land(j,k);
            elseif i ==4
                total(i,j,k) = f121994_land(j,k);
            elseif i ==5
                total(i,j,k) = f121995_land(j,k);
            elseif i ==6
                total(i,j,k) = f121996_land(j,k);
            elseif i ==7
                total(i,j,k) = f121997_land(j,k);
            elseif i ==8
                total(i,j,k) = f121998_land(j,k);
            elseif i ==9
                total(i,j,k) = f121999_land(j,k);
            elseif i ==10
                total(i,j,k) = f141997_land(j,k);
            elseif i ==11
                total(i,j,k) = f141998_land(j,k);
            elseif i ==12
                total(i,j,k) = f141999_land(j,k);
            elseif i ==13
                total(i,j,k) = f142000_land(j,k);
            elseif i ==14
                total(i,j,k) = f142001_land(j,k);
            elseif i ==15
                total(i,j,k) = f142002_land(j,k);
            elseif i ==16
                total(i,j,k) = f142003_land(j,k);
            elseif i ==17
                total(i,j,k) = f152000_land(j,k);
            elseif i ==18
                total(i,j,k) = f152001_land(j,k);
            elseif i ==19
                total(i,j,k) = f152002_land(j,k);
            elseif i ==20
                total(i,j,k) = f152003_land(j,k);
            elseif i ==21
                total(i,j,k) = f152004_land(j,k);
            elseif i ==22
                total(i,j,k) = f152005_land(j,k);
            elseif i ==23
                total(i,j,k) = f152006_land(j,k);
            elseif i ==24
                total(i,j,k) = f152007_land(j,k);
            elseif i ==25
                total(i,j,k) = f162004_land(j,k);
            elseif i ==26
                total(i,j,k) = f162005_land(j,k);
            elseif i ==27
                total(i,j,k) = f162006_land(j,k);
            elseif i ==28
                total(i,j,k) = f162007_land(j,k);
            elseif i ==29
                total(i,j,k) = f162008_land(j,k);
            elseif i ==30
                total(i,j,k) = f162009_land(j,k);
            elseif i ==31
                total(i,j,k) = f182010_land(j,k);
            end                
        end
    end
end

for k = 1:31
    open urban_EA.fig
    temp = zeros(size_i,size_j);
    for i = 1:size_i
        for j = 1:size_j
            temp(i,j) = (total(k,i,j)/63)*100;
        end
    end
    image(temp)
    hold on
    spy(data_mask,'k')
    hold on
end
clear i j


%% carbon map

total_C = zeros(31,size_i,size_j);
c_land = zeros(31,1);

for k = 1:31
    for i = 1:size_i
        for j = 1:size_j
            if (i>=1230)&&(j>1090)
            elseif (i>=1360)&&(j>865)
            elseif (i>775)&&((j>540)&&(j<1300))
                total_C(k,i,j) = (total(k,i,j)/63)*100;
            end
        end
    end
end

for k = 1:31
    n_land = 0;
    for i = 1:size_i
        for j = 1:size_j
                if total_C(k,i,j) >95
                    n_land = n_land + 1;
%                 else
%                     total_C(k,i,j) = 0;
                end
        end
    end
    c_land(k,1) = (carbon(k,1)/n_land)*10^12;
    for i = 1:size_i
        for j = 1:size_j
            if total_C(k,i,j)>95
                total_C(k,i,j) = c_land(k,1);
            elseif total_C(k,i,j)<95
                total_C(k,i,j) = 0;
            end
        end
    end     
end


%% piture for C map

temp = zeros(size_i, size_j);

for i = 1:31
    for j = 1:size_i
        for k = 1:size_j
            temp(j,k) = total_C(i,j,k);
        end
    end
    open carbon.fig
    image(temp)
    hold on
    spy(data_mask,'k')
    hold on
end
