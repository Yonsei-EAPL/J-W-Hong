%% History
% 160127 Mr. Je-Woo Hong
    % spectrum analysis with ts_data
    % including double rotation, mean wind-speed, wind-direction, and spectra analysis (u,v,w,Ts,CO2)
%% load input data
% dataDir = 'E:\EAPL\JW_Observation\ÀºÆò´ºÅ¸¿î\_csv'; % folder name
% dataName = dir(fullfile(dataDir, 'CSV_3302.ts_data_20*'))'; 

dataDir = 'e:\EAPL\JW_Observation\´öÀ¯»ê\20160322 È«Áø±Ô±³¼ö´Ô(´öÀ¯»êÇÃ·°½º_6¿ù)\20150716'; % folder name % EP
dataName = dir(fullfile(dataDir, 'CSV_5498.ts_data_*'))'; % Mt.´öÀ¯

% dataDir = 'h:\b_EAPL\JW_Observation\¿ÀÃ¢\_csv\'; % folder name
% dataName = dir(fullfile(dataDir, 'CSV_8665.ts_data_*'))'; %for ochang
% dataName = dir(fullfile(dataDir, 'CSV_7679.ts_data_*'))'; %for hongcheon
%dataName = dir(fullfile(dataDir, 'CSV_7681.ts_data_*'))'; %for samcheok
%dataName = dir(fullfile(dataDir, 'CSV_7682.ts_data_*'))'; %for pyeongchang
%dataName = dir(fullfile(dataDir, 'CSV_6330.ts_data_*'))'; %for jeju

total_result=cell(length(dataName),2);

%% position
po_u = 5;
po_v = 6;
po_w = 7;
po_Ts = 8;
po_CO2 = 10;
po_H2O = 11;
po_cell_prs = 14;
po_diff_prs = 17;

% % 2015-10 
% po_u = 6;
% po_v = 7;
% po_w = 8;
% po_Ts = 9;
% po_CO2 = 11;
% po_H2O = 12;
% po_cell_prs = 15;
% po_diff_prs = 18;



%% constant
R = 8.3143*10^(-6); % kPa*m3/K/mmol
Mc = 12; % mgC/mmol
Md = 0.029; % g/mmol
Mv = 0.018015; % g/mmol
% sonic_ang = 220+8.06; %for EP NewTown
sonic_ang = 0; %for Mt. ´öÀ¯
% sonic_ang = 230+8.04; %for SF
% sonic_ang = 90+7.24; %for BS Tower

%% main process
for nlm=1:length(dataName)
    data=importdata(fullfile(dataDir, dataName(nlm).name));    
    
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
%         if ((mod(data(i,4),100)==0)||(mod(data(i,4),100)==30))&&(mod(data(i,5),1)==0.1) % 2015-09 ~ 10-04
%         if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==0.05) %
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


    %% result for one-file (temp.)
    result = zeros(max(num_30min_n(:,1)),418);
    % 1; n_data (unitless)
    % 2; mean wind speed, U (m/s)
    % 3; wind direction (including sonic_angle) (degree)

    for i = 1:max(num_30min_n(:,1))
        result(i,1) = num_30min_n(i,2);    
        if num_30min_n(i,2)==18000
            % extract
            temp = zeros(18000,5);
            for j = 1:18000
                if i == 1
                    temp(j,1) = data(j,po_u);
                    temp(j,2) = data(j,po_v);
                    temp(j,3) = data(j,po_w);
                    temp(j,4) = data(j,po_Ts);
                    temp(j,5) = data(j,po_CO2);
                    temp(j,6) = data(j,po_H2O);
                    temp(j,7) = data(j,po_cell_prs);
                    temp(j,8) = data(j,po_diff_prs);
                else
                    temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
                    temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
                    temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
                    temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
                    temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
                    temp(j,6) = data(num_30min_n(i-1,3)+j,po_H2O);
                    temp(j,7) = data(num_30min_n(i-1,3)+j,po_cell_prs);
                    temp(j,8) = data(num_30min_n(i-1,3)+j,po_diff_prs);
                end
            end
            clear j

            % mean wind-speed, U
            temp_ws = 0;
            for j = 1:18000
                temp_ws = temp_ws + (temp(j,1)^2 + temp(j,2)^2 + temp(j,3)^2)^(0.5);
            end
            temp_ws = temp_ws/18000;
            result(i,2) = temp_ws; % result : mean wind-speed, U
            clear j %temp_ws

            % double rotation
            temp_wind = zeros(18000,3);
            u_bar = mean(temp(:,1));
            v_bar = mean(temp(:,2));
            w_bar = mean(temp(:,3));            
            alpha = atan2(v_bar,u_bar);
            beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.5)));
            for j = 1:18000
                temp_wind(j,1) = cos(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+sin(beta)*temp(j,3);
                temp_wind(j,2) = -sin(alpha)*temp(j,1)+cos(alpha)*temp(j,2);
                temp_wind(j,3) = -sin(beta)*(cos(alpha)*temp(j,1)+sin(alpha)*temp(j,2))+cos(beta)*temp(j,3);
                temp(j,1) = temp_wind(j,1);
                temp(j,2) = temp_wind(j,2);
                temp(j,3) = temp_wind(j,3);
            end
            clear u_bar v_bar w_bar alpha beta j temp_wind
            w_bar = mean(temp(:,3));     
            
            % Spectra analysis for u
            x = temp(:,1);
            nx=4;   % the number of variables
            ii=3;    %
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
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            % spectrum smoothing (high frequency band)
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
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+4) = fpsd_s(ii,1);
            end
            clear ii fpds_s

            % Spectra analysis for v
            x = temp(:,2);
            nx=4;   % the number of variables
            ii=3;    %
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
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            % spectrum smoothing (high frequency band)
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
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+50) = fpsd_s(ii,1);
            end
            clear ii fpds_s        

            % Spectra analysis for w
            x = temp(:,3);
            nx=4;   % the number of variables
            ii=3;    %
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
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            % spectrum smoothing (high frequency band)
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
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+96) = fpsd_s(ii,1);
            end
            clear ii fpds_s                

            % Spectra analysis for Ts
            x = temp(:,4);
            nx=4;   % the number of variables
            ii=3;    %
            dt=0.1;  % time difference between two conjecutive points (seconds)
            nst=16; % starting point of non-overlapping average
            pt_dec=10;  % how many points will be overlapped (point per decade)
            % Calculation of statistics
            mm=size(x);
            n=mm(1);
            n2=n/2;
            nvar=mm(2);
            xm(1,1)=mean(x);
            xstd(1,1)=std(x);
            xvar(1,1)=xstd(1,1).^2;
            xp(:,1)=x(:,1)-xm(1,1); % perturbation value
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            % spectrum smoothing (high frequency band)
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
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+142) = fpsd_s(ii,1);
            end
            clear ii fpds_s                

            % convert from mole fraction to density
%             Ts = mean(temp(:,4)); % result : mean_Ts 
%             Xc = mean(temp(:,5)); % result : mean co2
%             Xv = mean(temp(:,6)); % result : mean h2o
%             Tk = (Ts + 273.15)/(1+0.000321*Xv); % K ; AmeriFlux
%             amb_prs = (mean(temp(:,7))-mean(temp(:,8))); % kPa
%             rho_v = Xv*amb_prs*Mv/(R*Tk*(1000+Xv)); % g/m3
%             rho_c = Xc*Mc*10^(-6)*(amb_prs/(R*Tk) - rho_v/Mv); % mgC/m3
            for hjw= 1:18000
                Xv_temp = temp(hjw,6);
                Xc_temp = temp(hjw,5);
                Tk_temp = (temp(hjw,4) + 273.15)/(1+0.000321*Xv_temp); % K ; AmeriFlux
                amb_prs_temp = (temp(hjw,7)-temp(hjw,8)); % kPa
                temp(hjw,6) = Xv_temp * amb_prs_temp *Mv / (R*Tk_temp*(1000+Xv_temp)); % g/m3
                temp(hjw,5) = Xc_temp*Mc*10^(-6)*(amb_prs_temp/(R*Tk_temp) - temp(hjw,6)/Mv); % mgC/m3
            end
            clear hjw Xv_temp Xc_temp Tk_temp amb_prs_temp
            
            % Spectra analysis for CO2
            x = temp(:,5);
            nx=4;   % the number of variables
            ii=3;    %
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
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            % spectrum smoothing (high frequency band)
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
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+188) = fpsd_s(ii,1);
            end
            clear ii fpds_s                        
            
            % Spectra analysis for H2O
            x = temp(:,6);
            nx=4;   % the number of variables
            ii=3;    %
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
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            % spectrum smoothing (high frequency band)
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
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+234) = fpsd_s(ii,1);
            end
            clear ii fpds_s                                    
            
           %% cospectrum for w'T'
            fd2x = zeros(48,1);
            fd2y = zeros(48,1);           
            x = temp(:,3);
            y = temp(:,4);
            % Program History
            % NOV/30/2009   Original code (Jinkyu Hong)
            % Parameter setting
            nx=4;   % the number of variables
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
            ym(1,1)=mean(y);
            ystd(1,1)=std(y);
            yvar(1,1)=ystd(1,1).^2;
            yp(:,1)=y(:,1)-ym(1,1);
            coxy = 0;
            for j = 1:18000
                coxy = coxy + xp(j,1)*yp(j,1);
            end
            clear j
            coxy = coxy/18000;
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            fy(:,1)=fft(yp(:,1));
            fyr(:,1)=real(fy(:,1))./real(n);
            fyi(:,1)=imag(fy(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            co(:,1) = fxr(:,1).*fyr(:,1) + fxi(:,1).*fyi(:,1);
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            copsd = zeros(n2,1);
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
                copsd(j,1) = (2.*n.*dt)*co(j+1,1)./coxy;
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            copsd(n2,1) = (n.*dt)*co(n2+1,1)./coxy;
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            cofpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            cofpsd(:,1) = f(:).*copsd(:,1);    
            % spectrum smoothing (high frequency band)
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
            cofpsd_s=zeros(n2,1); % je-woo 2012-11-17    
            for js = 1:1:nst-1
                fpsd_s(js,1) = fpsd(js,1);
                cofpsd_s(js,1) = cofpsd(js,1);
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
                cosp = 0.0;
                while (flag1 == 'F')
                    p = p+1;
                    if(logx(p)<hi)&&(p<n2)
                        sp = sp + fpsd(p+1,1);
                        cosp = cosp + cofpsd(p+1,1);
                        fr = fr + f(p+1);
                        index = index + 1;
                        flah1 = 'F';
                    else
                        if(index>0)
                            k = k+1;
                            fpsd_s(k,1) = sp./real(index);
                            cofpsd_s(k,1) = cosp./real(index);
                            ax(k) = fr./real(index);
                        end
                        p = p-1;
                        flag1 = 'T';
                    end
                end
            end  
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            clear y ym yp ystd yvar fyr fyi fy cosp copsd cofpsd co
            cofd2y=zeros(48,1);
            for id3 = 1:48
                fd2x(id3,1) = ax(id3,1);
                fd2y(id3,1) = fpsd_s(id3,1);
                cofd2y(id3,1) = cofpsd_s(id3,1);
            end
            clear id3 ax fpsd_s cofpsd_s 
            result(i,1+279) = coxy;
            for jw = 1:46
                result(i,jw+280) = cofd2y(jw,1);
            end
            clear jw
            clear jewoo coxy
            
           %% cospectrum for w'CO2'
            fd2x = zeros(48,1);
            fd2y = zeros(48,1);           
            x = temp(:,3);
            y = temp(:,5);
            % Program History
            % NOV/30/2009   Original code (Jinkyu Hong)
            % Parameter setting
            nx=4;   % the number of variables
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
            ym(1,1)=mean(y);
            ystd(1,1)=std(y);
            yvar(1,1)=ystd(1,1).^2;
            yp(:,1)=y(:,1)-ym(1,1);
            coxy = 0;
            for j = 1:18000
                coxy = coxy + xp(j,1)*yp(j,1);
            end
            clear j
            coxy = coxy/18000;
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            fy(:,1)=fft(yp(:,1));
            fyr(:,1)=real(fy(:,1))./real(n);
            fyi(:,1)=imag(fy(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            co(:,1) = fxr(:,1).*fyr(:,1) + fxi(:,1).*fyi(:,1);
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            copsd = zeros(n2,1);
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
                copsd(j,1) = (2.*n.*dt)*co(j+1,1)./coxy;
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            copsd(n2,1) = (n.*dt)*co(n2+1,1)./coxy;
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            cofpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            cofpsd(:,1) = f(:).*copsd(:,1);    
            % spectrum smoothing (high frequency band)
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
            cofpsd_s=zeros(n2,1); % je-woo 2012-11-17    
            for js = 1:1:nst-1
                fpsd_s(js,1) = fpsd(js,1);
                cofpsd_s(js,1) = cofpsd(js,1);
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
                cosp = 0.0;
                while (flag1 == 'F')
                    p = p+1;
                    if(logx(p)<hi)&&(p<n2)
                        sp = sp + fpsd(p+1,1);
                        cosp = cosp + cofpsd(p+1,1);
                        fr = fr + f(p+1);
                        index = index + 1;
                        flah1 = 'F';
                    else
                        if(index>0)
                            k = k+1;
                            fpsd_s(k,1) = sp./real(index);
                            cofpsd_s(k,1) = cosp./real(index);
                            ax(k) = fr./real(index);
                        end
                        p = p-1;
                        flag1 = 'T';
                    end
                end
            end  
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            clear y ym yp ystd yvar fyr fyi fy cosp copsd cofpsd co
            cofd2y=zeros(48,1);
            for id3 = 1:48
                fd2x(id3,1) = ax(id3,1);
                fd2y(id3,1) = fpsd_s(id3,1);
                cofd2y(id3,1) = cofpsd_s(id3,1);
            end
            clear id3 ax fpsd_s cofpsd_s 
            result(i,1+325) = coxy;
            for jw = 1:46
                result(i,jw+326) = cofd2y(jw,1);
            end
            clear jw
            clear jewoo coxy            
            
           %% cospectrum for w'q'
            fd2x = zeros(48,1);
            fd2y = zeros(48,1);           
            x = temp(:,3);
            y = temp(:,6);
            % Program History
            % NOV/30/2009   Original code (Jinkyu Hong)
            % Parameter setting
            nx=4;   % the number of variables
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
            ym(1,1)=mean(y);
            ystd(1,1)=std(y);
            yvar(1,1)=ystd(1,1).^2;
            yp(:,1)=y(:,1)-ym(1,1);
            coxy = 0;
            for j = 1:18000
                coxy = coxy + xp(j,1)*yp(j,1);
            end
            clear j
            coxy = coxy/18000;
            % FFT
            fx(:,1)=fft(xp(:,1));
            fxr(:,1)=real(fx(:,1))./real(n);
            fxi(:,1)=imag(fx(:,1))./real(n);
            fy(:,1)=fft(yp(:,1));
            fyr(:,1)=real(fy(:,1))./real(n);
            fyi(:,1)=imag(fy(:,1))./real(n);
            % Power spectrum 
            s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
            co(:,1) = fxr(:,1).*fyr(:,1) + fxi(:,1).*fyi(:,1);
            % Power spectrum density
            psd = zeros(n2,1); % zero padding: power spectrum density
            copsd = zeros(n2,1);
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
                copsd(j,1) = (2.*n.*dt)*co(j+1,1)./coxy;
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            copsd(n2,1) = (n.*dt)*co(n2+1,1)./coxy;
            % frequency weighted power spectrum
            fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            fpsd(:,1) = f(:).*psd(:,1);
            cofpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
            cofpsd(:,1) = f(:).*copsd(:,1);    
            % spectrum smoothing (high frequency band)
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
            cofpsd_s=zeros(n2,1); % je-woo 2012-11-17    
            for js = 1:1:nst-1
                fpsd_s(js,1) = fpsd(js,1);
                cofpsd_s(js,1) = cofpsd(js,1);
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
                cosp = 0.0;
                while (flag1 == 'F')
                    p = p+1;
                    if(logx(p)<hi)&&(p<n2)
                        sp = sp + fpsd(p+1,1);
                        cosp = cosp + cofpsd(p+1,1);
                        fr = fr + f(p+1);
                        index = index + 1;
                        flah1 = 'F';
                    else
                        if(index>0)
                            k = k+1;
                            fpsd_s(k,1) = sp./real(index);
                            cofpsd_s(k,1) = cosp./real(index);
                            ax(k) = fr./real(index);
                        end
                        p = p-1;
                        flag1 = 'T';
                    end
                end
            end  
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            clear y ym yp ystd yvar fyr fyi fy cosp copsd cofpsd co
            cofd2y=zeros(48,1);
            for id3 = 1:48
                fd2x(id3,1) = ax(id3,1);
                fd2y(id3,1) = fpsd_s(id3,1);
                cofd2y(id3,1) = cofpsd_s(id3,1);
            end
            clear id3 ax fpsd_s cofpsd_s 
            result(i,1+371) = coxy;
            for jw = 1:46
                result(i,jw+372) = cofd2y(jw,1);
            end
            clear jw
            clear jewoo coxy
        end

    end
    total_result{nlm,1}= dataName(nlm).name;
    total_result{nlm,2}= result;
    clear data result
end
clear nlm