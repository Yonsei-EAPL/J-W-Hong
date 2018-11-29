%H:\b_MATLAB\2010_대학원\netcdf\20c3m\tas\ukmo_hadcm3\run2\

tas = nc_varget('tas_A1.nc','tas');
lat = nc_varget('tas_A1.nc','lat');
lon = nc_varget('tas_A1.nc','lon');

tas_n = zeros(1680,73,96);
for i = 1:73
    for j = 1:96
        temp1 = tas(:,i,j);
        temp2 = 0;
        
        for k = 1:12
            for l = 1:140
                temp2 = temp2 + temp1((l-1)*12+k,1);
            end
            temp2 = temp2/140;
            for l = 1:140
                temp1((l-1)*12+k,1) = temp1((l-1)*12+k,1) - temp2;
            end
        end
        clear k l
        
        for k = 1:1680
            tas_n(k,i,j) = temp1(k,1);
        end
        clear k
    end
end
clear i j

x=zeros(73,96);
for i = 1:73
    for j = 1:96
        x(i,j) = tas(20,i,j);
    end
end
clear i j
figure()
imagesc(x);figure(gcf);

x=zeros(73,96);
for i = 1:73
    for j = 1:96
        x(i,j) = tas_n(20,i,j);
    end
end
clear i j
figure()
imagesc(x);figure(gcf);



% for i = 1:73
%     for j = 1:96
%     x=tas_n(:,i,j);    
%     pause
%     wt(x)
%     end
% end
% clear i




%% 

data_set = zeros(1680,2);
data_set(:,1) = tas_n(:,52,35);
data_set(:,2) = tas_n(:,52,40);

%option 
num_data_day=12;  % year 
period_avg=140;
timestep=7200;
max_tau=24;
cycle=840;
% data_set information
size_data=size(data_set);
num_var=size_data(1,2);
num_data=size_data(1,1);
% clear size_data

max_bin = 600;
H=zeros(max_bin, num_var);
mi_result=zeros(max_bin, num_var, num_var);
for nbin=1:max_bin

    num_bin=nbin;

    % anomaly
    ano_data_set=anomaly(data_set, num_var, num_data, period_avg, timestep);
%     clear period_avg timestep num_data
    % anomaly information
    size_data=size(ano_data_set);
    num_ano_data=size_data(1,1);
    num_day=num_ano_data(1,1)/num_data_day;
    % clear size_data num_data_day

    % 130426_ anomaly cancel for electricity experiment
    % ano_data_set = data_set;

    % shannon entropycl
    for i=1:num_var
        H(nbin, i)=shannonentropy(ano_data_set, i, num_bin, num_ano_data);
    end

    % mutual information
    for i=1:num_var
        for j=1:num_var
            mi_result(nbin, i,j)=mutualinformation(ano_data_set, i, j, num_bin, num_ano_data);
        end
    end
 end

temp = zeros(max_bin,1);
temp2 = zeros(max_bin,1);
for i = 1:max_bin
    temp(i,1) = mi_result(i,1,1);
    temp2(i,1) = mi_result(i,2,1);
end
plot(temp)
hold on 
plot(temp2,'r')
hold on
for i = 1:max_bin
    temp2(i,1) = temp2(i,1)/temp(i,1);
end
plot(temp2,'c')




%%

x=tas_n(:,52,35);
%y=tas_n(:,35,69); % 동태평양 엘리뇨지역
% y=tas_n(:,61,26); % 러시아
% y=tas_n(:,55,73); % 미국
y=tas_n(:,52,57); % 태평양
%y=tas_n(:,52,25); % 중국


% set(gca,'YTickLabel',{'1','6','12','24','36','48','60','120','240','360','480'},...
%     'YTick',[0 2.584962501 3.584962501 4.584962501 5.169925001 5.584962501 5.906890596 6.906890596 7.906890596 8.491853096 8.906890596],...
%     'XTickLabel',{'1860','1870','1880','1890','1900','1910','1920','1930','1940','1950','1955','1960','1965','1970','1975','1980','1985','1990','1995','2000'},...
%     'XTick',[1 120 240 360 480 600 720 840 960 1080 1140 1200 1260 1320 1380 1440 1500 1560 1620 1680]);
% grid on

x=x(1081:1680,1);
y=y(1081:1680,1);

figure()
wtc(x,y)
set(gca,'YTickLabel',{'1','6','12','24','36','48','60','120','240','360','480'},...
    'YTick',[0 2.584962501 3.584962501 4.584962501 5.169925001 5.584962501 5.906890596 6.906890596 7.906890596 8.491853096 8.906890596],...
    'XTickLabel',{'1950','1955','1960','1965','1970','1975','1980','1985','1990','1995','2000'},...
    'XTick',[1 60 120 180 240 300 360 420 480 540 600]);
set(gcf,'color','w');
grid on


