[size_n size_var] = size(data);

result = zeros(size_n,1);

po_RH = 2;
po_T = 1;
po_P = 3;
unit_RH = 2; % 1:%/100, 2:%
unit_T = 1; % 1:degreeC, 2:K
unit_P = 2; % 1:kPa, 2:hPa

for i = 1:size_n
    % extract
    RH = data(i,po_RH);
    if unit_RH ==2
        RH = RH/100;
    end
    T = data(i,po_T);
    if unit_T ==2
        T = T - 273.15;
    end
    P = data(i,po_P);
    if unit_P ==2
        P = P/10;
    end

    % calculation
    Rv = 461.495;
    Rd = 287.05;    
    es = 6.1365*exp(17.502*T/(240.97+T))*0.1;
    e = es*RH;
    Pd = P - e;
    rho_v = (e*1000)/(Rv*(T+273.15));
    rho_d = (Pd*1000)/(Rd*(T+273.15));    
    RH = RH*100;
    Cpm = 1859 + 0.13*RH + (T)*(0.193+0.00569*RH) + (T)^2*(0.001+0.000005*RH);
    Cpd = 1005+(((T+273.15)-250.03)^2/3364);    
    rhoCp = Cpm*rho_v + Cpd*rho_d;
    
    % save result
    RH = RH/100;
    if (RH<0)||(RH>100)
        result(i,1) = -99999;
    elseif (T<-30)||(T>40)
        result(i,1) = -99999;
    elseif (P<90)||(P>120)
        result(i,1) = -99999;
    else
        result(i,1) = rhoCp;
    end
    
end
clear size_n size_var po_RH po_T po_P unit_RH unit_T unit_P
clear RH T P 
clear Rv Rd es e Pd rho_v rho_d Cpm Cpd rhoCp
clear i 