%% constant
R = 8.3143*10^(-6); % kPa*m3/K/mmol
Mc = 12; % mgC/mmol
Mco2 = 0.04401; % mg/umol
Md = 0.029; % g/mmol
Mv = 0.018015; % g/mmol
sonic_ang = 230+8.04; %for SF
% sonic_ang = 220+8.06; %for EP NewTown
% sonic_ang = 0; % Mr.´öÀ¯

%%
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
clear i num_30min

result = zeros(max(num_30min_n(:,1)),45*3+1);
result2 = zeros(max(num_30min_n(:,1)),20);
for ii = 1:max(num_30min_n(:,1))
    result(ii,1) = num_30min_n(ii,2);
    result2(ii,1) = num_30min_n(ii,2);    
    if num_30min_n(ii,2)==18000
        temp = zeros(num_30min_n(ii,2),5);

        for j = 1:num_30min_n(ii,2)
            if ii == 1
                temp(j,1) = data(j,8); % Ts
                temp(j,2) = data(j,26); % Xv 7200 [ mmol mol-1 ]
                temp(j,3) = data(j,37); % amb_prs 7500
                temp(j,4) = data(j,25); % CO2 7200 [ umol mol-1 ]
                temp(j,5) = data(j,33); % CO2 7500 [ mmol m-3 ]
                temp(j,6) = data(j,5); % u
                temp(j,7) = data(j,6); % v
                temp(j,8) = data(j,7); % w
                temp(j,9) = data(j,34); % H2O 7500 [ mmol m-3 ]
            else
                temp(j,1) = data(num_30min_n(ii-1,3)+j,8);
                temp(j,2) = data(num_30min_n(ii-1,3)+j,26);
                temp(j,3) = data(num_30min_n(ii-1,3)+j,37);
                temp(j,4) = data(num_30min_n(ii-1,3)+j,25);
                temp(j,5) = data(num_30min_n(ii-1,3)+j,33);
                temp(j,6) = data(num_30min_n(ii-1,3)+j,5); % u
                temp(j,7) = data(num_30min_n(ii-1,3)+j,6); % v
                temp(j,8) = data(num_30min_n(ii-1,3)+j,7); % w
                temp(j,9) = data(num_30min_n(ii-1,3)+j,34);                
            end
        end
        clear j

        % double rotation
        temp_wind = zeros(num_30min_n(ii,2),3);
        u_bar = mean(temp(:,6));
        v_bar = mean(temp(:,7));
        w_bar = mean(temp(:,8));
        alpha = atan2(v_bar,u_bar);
        beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.25)));
        for j = 1:num_30min_n(ii,2)
            temp_wind(j,1) = cos(beta)*(cos(alpha)*temp(j,6)+sin(alpha)*temp(j,7))+sin(beta)*temp(j,8);
            temp_wind(j,2) = -sin(alpha)*temp(j,6)+cos(alpha)*temp(j,7);
            temp_wind(j,3) = -sin(beta)*(cos(alpha)*temp(j,6)+sin(alpha)*temp(j,7))+cos(beta)*temp(j,8);
            temp(j,6) = temp_wind(j,1);
            temp(j,7) = temp_wind(j,2);
            temp(j,8) = temp_wind(j,3);
        end
        result2(ii,2)= alpha; % result : 1st rotation angle
        result2(ii,3)= beta; % result : 2nd rotation angle
        result2(ii,4)=std(temp(:,6)); % result : standard deviation of u
        result2(ii,5)=std(temp(:,7)); % result : standard deviation of v
        result2(ii,6)=std(temp(:,8)); % result : standard deviation of w
        result2(ii,7)=skewness(temp(:,6)); % result : skewness of u
        result2(ii,8)=kurtosis(temp(:,6)); % result : kurtosis of u
        result2(ii,9)=skewness(temp(:,7)); % result : skewness of v
        result2(ii,10)=kurtosis(temp(:,7)); % result : kurtosis of v
        result2(ii,11)=skewness(temp(:,8)); % result : skewness of w
        result2(ii,12)=kurtosis(temp(:,8)); % result : kurtosis of w    
        clear alpha beta j temp_wind        
        % mean wind-speed, U
        temp_ws = 0;
        for j = 1:num_30min_n(ii,2)
            temp_ws = temp_ws + (temp(j,6)^2 + temp(j,7)^2 + temp(j,8)^2)^(0.5);
        end
        temp_ws = temp_ws/num_30min_n(ii,2);
        result2(ii,13) = temp_ws; % result : mean wind-speed, U
        clear j temp_ws        
        % friction velocity
        uw = 0; % for mean(u'w')
        vw = 0; % for mean(v'w')
        TKE = 0; % for TKE
        for j = 1:num_30min_n(ii,2)
            uw = uw + (temp(j,6)-u_bar)*(temp(j,8)-w_bar);
            vw = vw + (temp(j,7)-v_bar)*(temp(j,8)-w_bar);
            TKE = TKE + (temp(j,6)-u_bar)^2+(temp(j,8)-w_bar)^2+(temp(j,7)-v_bar)^2;
        end
        uw = uw/num_30min_n(ii,2);
        vw = vw/num_30min_n(ii,2);
        TKE = TKE/2/num_30min_n(ii,2);
        result2(ii,14) = (uw^2 + vw^2)^(0.25); % result : u_star
        result2(ii,15) = uw; % result : cov_uw
        result2(ii,16) = uw/result2(ii,4)/result2(ii,6); % uw/stdu/stdw
        result2(ii,17) = TKE; % result : TKE
        clear uw vw TKE j        
        % std_u,v,w/u*,U
        result2(ii,18) = result2(ii,4)/result2(ii,14); % result : std_u/u*
        result2(ii,19) = result2(ii,5)/result2(ii,14); % result : std_v/u*
        result2(ii,20) = result2(ii,6)/result2(ii,14); % result : std_w/u*
        result2(ii,21) = result2(ii,4)/result2(ii,13); % result : std_u/U
        result2(ii,22) = result2(ii,5)/result2(ii,13); % result : std_v/U
        result2(ii,23) = result2(ii,6)/result2(ii,13); % result : std_w/U
       
        Ts = mean(temp(:,1));
        Xv = mean(temp(:,2));
        Xc = mean(temp(:,4));
        Tk = (Ts + 273.15)/(1+0.000321*Xv); % K ; AmeriFlux    
        result2(ii,24) = Tk-273.15; % result : actual temperature in degree C            
        amb_prs = mean(temp(:,3));
        V = 0.082*Tk*10^(-3)/(amb_prs/101.325); % m3/mol        
        e = Xv*amb_prs/(1000+Xv); % kPa
        es = 6.1365*exp(17.502*(Tk-273.15)/(240.97+(Tk-273.15)))/10; % kPa    
        RH = e/es*100;         
        rho_v = Xv*amb_prs*Mv/(R*Tk*(1000+Xv)); % g/m3        
        rho_c = Xc*Mc*10^(-6)*(amb_prs/(R*Tk) - rho_v/Mv); % mgC/m3        
        rho_co2 = Xc*10^(-6)*Mco2*(amb_prs/(R*Tk) - rho_v/Mv)*1000; % mg/m3        
        rho_d = (amb_prs-e)*Md/(R*Tk); % g/m3
        Cpd = 1005+((Tk-250.03)^2/3364); % J/K/kg ; AmeriFlux
        Cpm = 1859 + 0.13*RH + (Tk-273.15)*(0.193+0.00569*RH) + (Tk-273.15)^2*(0.001+0.000005*RH); % J/K/kg ; AmeriFlux
        rhoCp = Cpm*(rho_v/1000) + Cpd*(rho_d/1000); % specific heat capacity (J/K/kg)        
        cov_wTs = 0;
        cov_wq7200 = 0;
        cov_wq7500 = 0;
        V_bar = mean(temp(:,9));        
        for j = 1:num_30min_n(ii,2)
            cov_wTs = cov_wTs + (temp(j,8)-w_bar)*(temp(j,1)-Ts);
            cov_wq7200 = cov_wq7200 + (temp(j,8)-w_bar)*(temp(j,2)-Xv); % mmol/mol *m/s
            cov_wq7500 = cov_wq7500 + (temp(j,8)-w_bar)*(temp(j,9)-V_bar); % mmol/m3 * m/s
        end
        clear j        
        cov_wTs = cov_wTs/num_30min_n(ii,2);
        cov_wq7200 = cov_wq7200/num_30min_n(ii,2);
        cov_wq7500 = cov_wq7500/num_30min_n(ii,2);
        Qh = rhoCp * (cov_wTs - (0.000321*Tk*cov_wq7200)); % AmeriFlux
        result2(ii,25) = Qh;
        lambda = 2500.8 - 2.3668*Ts; % J/g
        V = 0.082*Tk*10^(-3)/(amb_prs/101.325); % m3/mol
        Qe7200 = lambda*Mv*cov_wq7200/V; 
        result2(ii,26) = Qe7200;
        Qe7500 = lambda*Mv*cov_wq7500;
        result2(ii,27) = Qe7500;
        co2_7200 = mean(temp(:,4));
        co2_7500 = mean(temp(:,5));
        cov_wc7200 = 0;
        cov_wc7500 = 0;
        for j = 1:num_30min_n(ii,2)
            cov_wc7200 = cov_wc7200 + (temp(j,8)-w_bar)*(temp(j,4)-co2_7200);
            cov_wc7500 = cov_wc7500 + (temp(j,8)-w_bar)*(temp(j,5)-co2_7500);
        end
        clear j 
        cov_wc7200 = cov_wc7200/num_30min_n(ii,2); % umol/mol * m/s
        cov_wc7500 = cov_wc7500/num_30min_n(ii,2); % mmol/m3 * m/s
        cov_wc7200 = cov_wc7200/V; % umol/m2/s 
        result2(ii,28) = cov_wc7200;
        result2(ii,29) = cov_wc7500*1000;
        % WPL corection for latent heat flux
        u = 1.6077; % ratio of molar masses of air to water
        cov_wq7500 = cov_wq7500*Mv; % g/m2/s
        E = (1+u*(rho_v/rho_d))*(cov_wq7500 + (Qh/rhoCp) * (rho_v/Tk)); % g/m2/s
        Qe_WPL = lambda*E; % W/m2
        result2(ii,30) = Qe_WPL;
        % WPL corection for co2 flux
        Fc_WPL = cov_wc7500*1000*Mco2 + u*(E/rho_d)*(rho_co2/(1+u*(rho_v/rho_d))) + ((Qh/rhoCp)*(rho_co2)/Tk); % mg/m2/s
        Fc_WPL = Fc_WPL/Mco2; % umol/m2/s
        result2(ii,31) = Fc_WPL;
        clear u_bar v_bar w_bar Qh Qe7200 Qe7500 co2_7200 co2_7500 lambda V_bar cov_wTs cov_wq7200 cov_wq7500 cov_wc7200 cov_wc7500
        clear Ts Xv Xc Tk amb_prs rho_v rho_c rho_co2 rho_d Cpd Cpm rhoCp e es RH u E Qe_WPL Fc_WPL
        
        % spectrum
        for j = 1:num_30min_n(ii,2)
            temp(j,5) = temp(j,5)*1000;
            temp(j,4) = temp(j,4)/V;
        end
        clear j V
        
        x = temp(:,4);
        nx=4;   % the number of variables
        i=3;    %
        dt=0.1;  % time difference between two conjecutive points (seconds)
        nst=16; % starting point of non-overlapping average
        pt_dec=10;  % how many points will be overlapped (point per decade)
        mm=size(x);
        n=mm(1);
        n2=n/2;
        nvar=mm(2);
        xm(1,1)=mean(x);
        xstd(1,1)=std(x);
        xvar(1,1)=xstd(1,1).^2;
        xp(:,1)=x(:,1)-xm(1,1); % perturbation value
        fx(:,1)=fft(xp(:,1));
        fxr(:,1)=real(fx(:,1))./real(n);
        fxi(:,1)=imag(fx(:,1))./real(n);
        s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
        psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
        f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
        for j=1:1:n2-1
            f(j) = real(j)/(real(n).*dt);
            psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
        end
        f(n2) = real(n2)/(real(n).*dt);
        psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
        display ('Test the Parseval Theorem')
        display ('The value below should be near around 1.0 for the energy conservation')
        sum(psd)./(n.*dt)
        fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
        fpsd(:,1) = f(:).*psd(:,1);
        logx=zeros(n2,1); % je-woo 2012-11-17
        for j = nst:1:n2
            logx(j) = log10(f(j));
        end
        lf_mi = logx(nst);
        lf_ma = logx(n2);
        nn = int8((lf_ma - lf_mi).*real(pt_dec));
        ax=zeros(n2,1); % je-woo 2012-11-17
        for js = 1:1:nst-1
            ax(js) = f(js);
        end
        fpsd_s=zeros(n2,1); % je-woo 2012-11-17
        for js = 1:1:nst-1
            fpsd_s(js,1) = fpsd(js,1);
        end
        p = nst - 1;
        k = nst - 1;
        sp = 0.0;
        for js = 1:1:45
            flag1 = 'F';
            lo = (lf_mi) + real(js-1)./real(pt_dec);
            hi = (lf_mi) + real(js)./real(pt_dec);
            fr = 0.0;
            index = 0;
            sp = 0.0;
            while (flag1 == 'F')
                p = p+1;
                if(logx(p)<hi)&&(p<n2)
                    sp = sp + fpsd(p+1,1);
                    fr = fr + f(p+1);
                    index = index + 1;
                    flah1 = 'F';
                else
                    if(index>0)
                        k = k+1;
                        fpsd_s(k,1) = sp./real(index);
                        ax(k) = fr./real(index);
                    end
                    p = p-1;
                    flag1 = 'T';
                end
            end
        end  
        clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
        clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
        fpsd_s7200 = fpsd_s;
        ax7200 = ax;
        
        x= temp(:,5);
        nx=4;   % the number of variables
        i=3;    %
        dt=0.1;  % time difference between two conjecutive points (seconds)
        nst=16; % starting point of non-overlapping average
        pt_dec=10;  % how many points will be overlapped (point per decade)
        mm=size(x);
        n=mm(1);
        n2=n/2;
        nvar=mm(2);
        xm(1,1)=mean(x);
        xstd(1,1)=std(x);
        xvar(1,1)=xstd(1,1).^2;
        xp(:,1)=x(:,1)-xm(1,1); % perturbation value
        fx(:,1)=fft(xp(:,1));
        fxr(:,1)=real(fx(:,1))./real(n);
        fxi(:,1)=imag(fx(:,1))./real(n);
        s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
        psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
        f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
        for j=1:1:n2-1
            f(j) = real(j)/(real(n).*dt);
            psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
        end
        f(n2) = real(n2)/(real(n).*dt);
        psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
        display ('Test the Parseval Theorem')
        display ('The value below should be near around 1.0 for the energy conservation')
        sum(psd)./(n.*dt)
        fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
        fpsd(:,1) = f(:).*psd(:,1);
        logx=zeros(n2,1); % je-woo 2012-11-17
        for j = nst:1:n2
            logx(j) = log10(f(j));
        end
        lf_mi = logx(nst);
        lf_ma = logx(n2);
        nn = int8((lf_ma - lf_mi).*real(pt_dec));
        ax=zeros(n2,1); % je-woo 2012-11-17
        for js = 1:1:nst-1
            ax(js) = f(js);
        end
        fpsd_s=zeros(n2,1); % je-woo 2012-11-17
        for js = 1:1:nst-1
            fpsd_s(js,1) = fpsd(js,1);
        end
        p = nst - 1;
        k = nst - 1;
        sp = 0.0;
        for js = 1:1:45
            flag1 = 'F';
            lo = (lf_mi) + real(js-1)./real(pt_dec);
            hi = (lf_mi) + real(js)./real(pt_dec);
            fr = 0.0;
            index = 0;
            sp = 0.0;
            while (flag1 == 'F')
                p = p+1;
                if(logx(p)<hi)&&(p<n2)
                    sp = sp + fpsd(p+1,1);
                    fr = fr + f(p+1);
                    index = index + 1;
                    flah1 = 'F';
                else
                    if(index>0)
                        k = k+1;
                        fpsd_s(k,1) = sp./real(index);
                        ax(k) = fr./real(index);
                    end
                    p = p-1;
                    flag1 = 'T';
                end
            end
        end  
        loglog(f(:),fpsd(:,1),'xk',ax(:),fpsd_s(:,1),'-or')
        clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
        clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
        y=[0.0100000000000000,1;10,0.0100000000000000;];
        hold on
        plot(y(:,1),y(:,2))
        fpsd_s7500 = fpsd_s;
        ax7500 = ax;

        for j= 1:45
            result(ii,1+j) = ax(j,1);
            result(ii,46+j) = fpsd_s7200(j,1);
            result(ii,91+j) = fpsd_s7500(j,1);
        end
        clear j
        
    end   
end
clear ii Mc Mco2 Md Mv R ax ax7200 ax7500 fpsd_s fpsd_s7200 fpsd_s7500 num_30min_n size_n size_var sonic_ang temp x y 

% loglog(ax7200,fpsd_s7200,'ko-')
% hold on
% loglog(ax7500,fpsd_s7500,'rx-')
% y=[0.0100000000000000,1;10,0.0100000000000000;];
% hold on
% plot(y(:,1),y(:,2))