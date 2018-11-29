%% Information
% To estimate fi_m (dimensionless wind shear) and
%             fi_h (dimensionless potential-temperature gradient)
%             using 3m & 10m CSAT data from EP in 2010, 2011, and 2012.
% Macdonald et al. and Li-DAR DEM are used for z0 and zd by 45 degree WD.
%             (Je-Woo Hong, 2014, SNU master course thesis)


%% history
% 2014-04-11 J-W Hong, 


%% input
% 1 : u*, friction velocity, (u'w'^2 + v'w'^2)^(0.5)
% 2 : Ts_10m, sonic temperature @ 10m (degreeC)
% 3 : Ts_3m,  sonic temperature @ 3m (degreeC)
% 4 : Q_h_10m, sensible heat flux @ 10m (Wm-2) 
% 5 : U_10m, wind-speed @ 10m (ms-1)
% 6 : U_3m, wind-speed @ 3m (ms-1)  
% 7 : WD_10m, wind-direction @ 10m (degree)


%% constant
zd_WD = [5.7;6.7;7.5;13.5;10.8;10.9;11.9;3.4]; % zero-plane displacement
z0_WD = [0.4;0.6;1.7;0.5;0.5;0.3;1.0;0.4]; % roughness length
k = 0.4; % von Karmann constant
g = 9.8; % gravity acceleration
zm = 30; % measurement height
dzm = 7;
rhoCp = 1230; % used from EDIRE program


%% main
temp = size(input);
size = temp(1,1);
clear temp 

result = zeros(size, 6);
% 1 : L, Obukhov length
% 2 : zoL
% 3 : fi_m, dimensionless wind shear
% 4 : fi_h, dimensionless potential-temperature gradient
% 5 : psi_m, integrated Businger-Dyer relationship for momentum 
% 6 : psi_h,

z_WD = zeros(8,1); % z = zm-zd
for i = 1:8
    z_WD(i,1) = zm - zd_WD(i,1);
end
clear i

for i = 1:size
    input(i,1) = input(i,1)^(0.5);
end
clear i

for i = 1:size
    L = -1*(input(i,1)^3)/(k*(g/(input(i,2)+273.15))*(input(i,4)/rhoCp));
    result(i,1) = L;
    
    WD_index = (input(i,7) - mod(input(i,7),45))/45 + 1;
  
    zoL = z_WD(WD_index,1)/L;
    result(i,2) = zoL;
    
    fi_m = (k*z_WD(WD_index,1)/input(i,1))*((input(i,5)-input(i,6))/dzm); 
    result(i,3) = fi_m;
    
    ts = -1*input(i,4)/(rhoCp*input(i,1)); % temperature scale
    fi_h = (k*z_WD(WD_index,1)/ts)*((input(i,2)-input(i,3))/dzm);
    result(i,4) = fi_h;
    
    psi_m = log(zm/z0_WD(WD_index,1)) - input(i,5)*0.4/input(i,1);
    result(i,5) = psi_m;
    
    psi_h = log(zm/z0_WD(WD_index,1)) - k/ts*(input(i,2)- ((input(i,2)-input(i,3))/dzm*(z0_WD(WD_index,1)-zm)));
    result(i,6) = psi_h;
end
clear i


%% clear 
clear k g rhoCp dzm zm zd_WD z0_WD
clear z_WD L WD_index zoL fi_m ts fi_h psi_m psi_h


%% plot
x = [-2.5:0.01:0.5]; % x(1,301)
y = zeros(301,1);
for i = 1:301
    if x(1,i) >=0
        y(i,1) = 0.95 + 7.8*x(1,i);
    else
        y(i,1) = 0.95*(1-11.6*x(1,i))^(-0.5);
    end
end
clear i
figure(1)
plot(x(1,:),y(:,1))
hold on
clear x y

temp = 0;
for i = 1:size
    if (input(i,5)>1.5)&&((input(i,5)>input(i,6))&&((result(i,2)>-2.5)&&(result(i,2)<0.5)))
        temp = temp + 1;
    end
end
temp_r = zeros(temp,2);
for i = 1:temp
    if (input(i,5)>1.5)&&((input(i,5)>input(i,6))&&((result(i,2)>-2.5)&&(result(i,2)<0.5)))
        temp_r(i,1) = result(i,2);
        temp_r(i,2) = result(i,3);
    end
end
clear i temp
figure(1)
plot(temp_r(:,1),temp_r(:,2),'xr')
clear temp_r
