%% history ------------------------------------------------------
% This program made by Je-Woo Hong (jewoo@yonsei.kr), in feb-2014
% (Ecosystem-Atmosphere Processes Lab., Yonsei University, South Korea)
%
% Object of this program is to estimate the surface flux using surface
% renewal method (following Castellvi et al (2008) A&FM)
%
% For this program, CSV type data from closed-path EC system is used.
% ------------------------------------------------------------------------
%% constant -----------------------------------------------------
zm = 30; % measurement height in meter
h = 13; % mean building height (or canopy height) in meter
    d = 0.7*h; % zero plane displacement (Brutsaert, 1988)
    z0 = 0.12*h; % roughness length (Brutsaert, 1988)
z = zm - d;
da = 8.06; % declination angle in degree (ex; EP 8.06/ SF 8.04)
g = 9.7991; % gravity acceleration (ex; EP 9.7991/ SF 9.7992) (KRISS)
k = 0.4; % von Karman constant
% ------------------------------------------------------------------------
%% input file (csv file) ----------------------------------------
% sonic_ang = 230+d; %for SF
sonic_ang = 220+da; %for EP NewTown
[size_n size_var] = size(data);
num_30min = 0;
num_30min_n = zeros(1,3);
for i = 1:size_n
    if i==1
        num_30min = 1;
        num_30min_n(num_30min,1) = num_30min;
        num_30min_n(num_30min,2) = 1;
        num_30min_n(num_30min,3) = 1;
    else
        if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==0.1)
            num_30min = num_30min+1;
            num_30min_n(num_30min,1) = num_30min;
            num_30min_n(num_30min,2) = 1;
            num_30min_n(num_30min,3) = num_30min_n(num_30min-1,3)+1;
        else
            num_30min_n(num_30min,2) = num_30min_n(num_30min,2)+1;
            num_30min_n(num_30min,3) = num_30min_n(num_30min,3)+1;
        end
    end
end
clear i 
% ------------------------------------------------------------------------
%% result -------------------------------------------------------
result = zeros(num_30min,12);
% 1; ur, mean horizontal wind speed
% 2; rx, time-lag (in second)
% 3; tau, time scale in seconds
% 4; A, amplitude in degree C
% 5; Ts, mean sonic(virtual) temperature 
% 6; u*, friction velocity from iteration
% 7; alpha, scale factor in SR method, from iteration
% 8; kinematic H, from iteration
% 9; L, obukhov length from iteration
% 10; zeta, zoL from iteration
% 11; temp_n, number of iteration (100 is maximum)
% 12; H, sensible heat flux from SR method
% ------------------------------------------------------------------------
%% position -----------------------------------------------------
po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8;
po_CO2 = 10;
po_H2O = 11;
po_cell_tmpr = 13;
po_cell_prs = 14;
% ------------------------------------------------------------------------
%% calculation --------------------------------------------------
for i = 1:max(num_30min_n(:,1))
    % extract ------------------------------------------------------------
    temp = zeros(num_30min_n(i,2),8);
    temp_T = zeros(num_30min_n(i,2),1);
    for j = 1:num_30min_n(i,2)
        if i ==1
            temp(j,1) = data(j,po_u);
            temp(j,2) = data(j,po_v);
            temp(j,3) = data(j,po_w);
            temp(j,4) = data(j,po_Ts);
            temp_T(j,1) = data(j,po_Ts);
            temp(j,5) = data(j,po_CO2);
            temp(j,6) = data(j,po_H2O);
            temp(j,7) = data(j,po_cell_tmpr);
            temp(j,8) = data(j,po_cell_prs);
        else
            temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
            temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
            temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
            temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
            temp_T(j,1) = data(num_30min_n(i-1,3)+j,po_Ts);
            temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
            temp(j,6) = data(num_30min_n(i-1,3)+j,po_H2O);
            temp(j,7) = data(num_30min_n(i-1,3)+j,po_cell_tmpr);
            temp(j,8) = data(num_30min_n(i-1,3)+j,po_cell_prs);
        end
    end
    clear j
    % --------------------------------------------------------------------
    % rx (time-lag in seconds) -------------------------------------------
    max_n = 20;
    temp_r = zeros(max_n,1);
    for j = 1:max_n
        temp_r(j,1) = abs(Sf(temp_T,3,j)/(0.1*j)); 
    end
    clear j
    [temp_r1 temp_r2] = max(temp_r);
    clear temp_r1
    rx = temp_r2; % time-lag in number
    clear temp_r2
    result(i,2) = rx*0.1;
    % --------------------------------------------------------------------
    % tau ----------------------------------------------------------------
    p = 10*Sf(temp_T,2,rx)-Sf(temp_T,5,rx)/Sf(temp_T,3,rx);
    q = 10*Sf(temp_T,3,rx);
    A = ((-1)*q/2 + ((q/2)^(2)+(p/3)^(3))^(0.5))^(1/3);
    A = A + ((-1)*q/2 - ((q/2)^(2)+(p/3)^(3))^(0.5))^(1/3);
    A = real(A);
    tau = (-1)*((A^(3))*(rx*0.1)/Sf(temp_T,3,rx));
    if tau<0
        tau = tau*(-1);
    end
    result(i,3) = tau;
    result(i,4) = A;
    clear p q temp_n
    % --------------------------------------------------------------------
    % mean horizontal wind speed -----------------------------------------
    temp_u = 0;
    for j = 1:num_30min_n(i,2)
        temp_u = temp_u + (temp(j,1)^2 + temp(j,2)^2)^(0.5);
    end
    ur = temp_u/num_30min_n(i,2);
    clear j
    clear temp_u
    result(i,1) = ur;
    % --------------------------------------------------------------------
    % mean Ts (virtual temperature) --------------------------------------
    Ts = mean(temp_T(:,1));
    result(i,5) = Ts;
    % --------------------------------------------------------------------
    % iteration ----------------------------------------------------------
    %   step1. start from neutral condition (zeta ; zoL = 0)
    %   step2. u*
    %   step3. alpha, H
    %   step4. L
    %   step5. zeta (first result) !
    %   threshold : delta u* < 0.005 m s-1
    zeta = 0;
    du_star = 1;
    temp_n = 0;
    while (du_star>0.0005)&&(temp_n < 500)
        temp_n = temp_n + 1;
        u_star = k*ur/log(z/z0) - 0.1*Psi(zeta);
        if u_star<0
            u_star = 0.01;
        end
        if temp_n == 1
            x_u_star = u_star;
        else
            du_star = abs(x_u_star - u_star);
            x_u_star = u_star;
        end
        alpha = (k/pi() * z/(zm^2) * tau * u_star / Fi(zeta))^(0.5);
        kH = alpha * zm * A / tau;
        L = (u_star)^3 / (k*g*kH) * (Ts+273.15);
        zeta = z/L;
    end
    result(i,6) = u_star;
    result(i,7) = alpha;
    result(i,8) = kH;
    result(i,9) = L;
    result(i,10) = zeta;
    result(i,11) = temp_n;
    % --------------------------------------------------------------------
    % H from SR method ---------------------------------------------------
    H = alpha * z * A / tau;
    result(i,12) = H;
    % --------------------------------------------------------------------
end
clear i temp
% ------------------------------------------------------------------------
%% clear --------------------------------------------------------
clear zm h d z0 z da g k sr % constant
clear sonic_ang size_n size_var num_30min num_30min_n % input file
clear po_u po_v po_w po_Ts po_CO2 po_H2O po_cell_tmpr po_cell_prs % position
clear L Ts  du_star kH max_n rx temp_T temp_n temp_r u_star ur x_u_star zeta % calculation
clear alpha A tau H % calculation
% ------------------------------------------------------------------------
