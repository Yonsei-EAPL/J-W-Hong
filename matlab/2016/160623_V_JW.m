%%

%% Velocity (V(D))
% x=[1.0e-06;2.0e-06;3.0e-06;4.0e-06;5.0e-06;6.0e-06;7.0e-06;8.0e-06;9.0e-06;1.0e-05;2.0e-05;3.0e-05;4.0e-05;5.0e-05;6.0e-05;7.0e-05;8.0e-05;9.0e-05;0.00010;0.00020;0.00030;0.00040;0.00050;0.00060;0.00070;0.00080;0.00090;0.0010;0.0020;0.0030;0.0040;0.0050;0.0060;0.0070;];
[length a] = size(x);
clear a
y = zeros(length,2);

% constant
rho = 1.204; % kg m-3 @ 20 oC at 1013 hPa
rho_i = 998.23; % kg m-3
drho = rho_i - rho;
l0 = 6.62* 10^(-6); % cm
    l0 = l0*0.01; % m
sigma = 0.07275; % N m-1 = kg s-2
eta = 0.0001818; % g cm-1
    eta = eta*0.001*100; % kg m-1
g = 9.8; % m s-2

b0 = -3.18657; % for regime2
b1 = 0.992696;
b2 = -0.153193*10^(-2); 
b3 = -0.987059*10^(-3);
b4 = -0.578878*10^(-3);
b5 = 0.855176*10^(-4);
b6 = -0.327815*10^(-5);

c0 = -5.00015; % for regime3
c1 = 5.23778;
c2 = -2.04914;
c3 = 0.475294;
c4 = -0.542819*10^(-1);
c5 = 0.238449*10^(-2);

for i = 1:length

    d0 = x(i,1);
    
    if d0< 2*10^(-5)
        % regime1
        C1 = drho*g/(18*eta);
        Csc = 1+2.51*(l0/d0);
        V = C1*Csc*d0^2;

    elseif d0< 1*10^(-3)
        % regime2
        Csc = 1+2.51*(l0/d0);
        C2 = 4*rho*drho*g/(3*eta^2);
        NDa = C2*d0^3;
        X = log(NDa);
        Y = b0+ b1*X+ b2*X^2+ b3*X^3+ b4*X^4+ b5*X^5+ b6*X^6;
        NRe = Csc*exp(Y);
        V = eta*NRe/(rho*d0);

    else
        % regime3
        Np = sigma^3*rho^2/(eta^4*drho*g);
        C3 = 4*drho*g/3/sigma;
        Bo = C3*d0^2;
        X = log(Bo*Np^(1/6));
        Y = c0+c1*X+c2*X^2+c3*X^3+c4*X^4+c5*X^5;
        NRe = Np^(1/6)*exp(Y);
        V = eta*NRe/(rho*d0);
    end
    
    y(i,1) = V;
end
clear i %length
clear rho drho l0 sigma eta g b0 b1 b2 b3 b4 b5 b6 c0 c1 c2 c3 c4 c5 d0 C1 Csc V C2 NDa X Y NRe Np C3 Bo %rho_i
% V=y;
% clear y;

% figure1 = figure('Color',[1 1 1]);
% axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
%     'XScale','log',...
%     'XMinorTick','on',...
%     'FontWeight','bold',...
%     'FontSize',20);
% box(axes1,'on');
% hold(axes1,'all');
% loglog(x,y(:,1),'MarkerSize',20,'Marker','*','LineWidth',3);
% xlabel('D  (m)','FontWeight','bold','FontSize',24);
% ylabel('V  (m s^{-1})','FontWeight','bold','FontSize',24);
% %ylabel('EXP  (inside of integral)','FontWeight','bold','FontSize',24);
% clear figure1 axes1

%% exponent

for i = 1:19
    for j = 1:63
        for k = 1:length
            E = x(k,1)^3*N0(j,1)*exp(-lambda(i,1)*x(k,1))*y(k,1);
            y(k,2)=E;
        end
        % polyarea
        x2 = zeros(length+2,1);
        y2 = zeros(length+2,1);

        x2(1,1) = x(1,1);
        y2(1,1) = 0;
        x2(length+2,1) = x(length,1);
        y2(length+2,1) = 0;
        for k = 1:length
            x2(1+k,1) = x(k,1);
            y2(1+k,1) = y(k,2);
        end
        area = polyarea(x2,y2);
        area = area*rho_i*pi()/6;
        y3(i,j) = area;
    end
end
clear i j k
clear E area rho_i 
clear length x2 y2 y


