%% History
% 160127 Mr. Je-Woo Hong
    % spectrum analysis with ts_data
    % including double rotation, mean wind-speed, wind-direction, and spectra analysis (u,v,w,Ts,CO2)
%% load input data
dataDir = 'E:\EAPL\JW_Observation\은평뉴타운\_csv'; % folder name
dataName = dir(fullfile(dataDir, 'CSV_3302.ts_data_2015_10_04*'))'; 
% dataDir = 'h:\b_EAPL\JW_Observation\오창\_csv\'; % folder name
% dataName = dir(fullfile(dataDir, 'CSV_8665.ts_data_*'))'; %for ochang
% dataName = dir(fullfile(dataDir, 'CSV_7679.ts_data_*'))'; %for hongcheon
%dataName = dir(fullfile(dataDir, 'CSV_7681.ts_data_*'))'; %for samcheok
%dataName = dir(fullfile(dataDir, 'CSV_7682.ts_data_*'))'; %for pyeongchang
%dataName = dir(fullfile(dataDir, 'CSV_6330.ts_data_*'))'; %for jeju

total_result=cell(length(dataName),2);

for nlm=1:length(dataName)
    data=importdata(fullfile(dataDir, dataName(nlm).name));    
    
%% input information
% sonic_ang = 230+8.04; %for SF
% sonic_ang = 90+7.24; %for BS Tower
sonic_ang = 220+8.06; %for EP NewTown
% sonic_ang =0;
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
        if ((mod(data(i,4),100)==0)||(mod(data(i,4),100)==30))&&(mod(data(i,5),1)==0.1)
%         if ((mod(data(i,3),100)==0)||(mod(data(i,3),100)==30))&&(mod(data(i,4),1)==0.05)
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


%% result
result = zeros(max(num_30min_n(:,1)),234+46);
% 1; n_data (unitless)
% 2; mean wind speed, U (m/s)
% 3; wind direction (including sonic_angle) (degree)


%% position
po_u = 6;
po_v = 7;
po_w = 8;
po_Ts = 9;
% po_CO2 = 11; %H2O분석에 사용
po_CO2 = 11;
% po_u = 7;
% po_v = 8;
% po_w = 9;
% po_Ts = 10;
% po_CO2 = 12;


%% main process
% wait = waitbar(0,'...Flux Data Process...');
jw_co = zeros(48,48);
% for jewoo = 49:96
    
    for i = 1:max(num_30min_n(:,1))
%         waitbar((i/(max(num_30min_n(:,1)))),wait,sprintf('%3f',i))

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
                else
                    temp(j,1) = data(num_30min_n(i-1,3)+j,po_u);
                    temp(j,2) = data(num_30min_n(i-1,3)+j,po_v);
                    temp(j,3) = data(num_30min_n(i-1,3)+j,po_w);
                    temp(j,4) = data(num_30min_n(i-1,3)+j,po_Ts);
                    temp(j,5) = data(num_30min_n(i-1,3)+j,po_CO2);
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
            clear j temp_ws

            % mean wind-direction
            temp_std_wd = 0;
            u_bar = mean(temp(:,1));
            v_bar = mean(temp(:,2));
            w_bar = mean(temp(:,3));
            if u_bar>0
                if v_bar>0
                    temp_wd = 360 - atan(v_bar/u_bar)/pi()*180;
                else
                    temp_wd = atan((-1*v_bar)/u_bar)/pi()*180;
                end
            else
                if v_bar>0
                    temp_wd = 180 + atan(v_bar/(-1*u_bar))/pi()*180;
                else
                    temp_wd = 180 - atan(v_bar/u_bar)/pi()*180;
                end
            end
            mean_wd = temp_wd;
            if mean_wd+sonic_ang>360
                mean_wd = mean_wd + sonic_ang-360;
            else
                mean_wd = mean_wd + sonic_ang;
            end
            result(i,3) = mean_wd; % result : wind direction
            clear j temp_wd temp_std_wd mean_wd 

            % double rotation
            temp_wind = zeros(18000,3);
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

            % Spectra analysis
            x = temp(:,1);
            nx=4;   % the number of variables
            ii=3;    %
            dt=0.1;  % time difference between two conjecutive points (seconds)
            nst=16; % starting point of non-overlapping average
            pt_dec=10;  % how many points will be overlapped (point per decade)
            % Calculation of statistics
            %t=[0.0:dt:1800-dt];
            %x(:,1)=xx(:,i);
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
            psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % Check the Parseval theorem
%             display ('Test the Parseval Theorem')
%             display ('The value below should be near around 1.0 for the energy conservation')
%             sum(psd)./(n.*dt)
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
            %for js = 1:1:nn
    %         for js = 1:1:28
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
            % Plotting
            %loglog(f(:),fpsd(:,1),'xk',ax(:),fpsd_s(:,1),'-or')
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+4) = fpsd_s(ii,1);
            end
            clear ii fpds_s

            % Spectra analysis
            x = temp(:,2);
            nx=4;   % the number of variables
            ii=3;    %
            dt=0.1;  % time difference between two conjecutive points (seconds)
            nst=16; % starting point of non-overlapping average
            pt_dec=10;  % how many points will be overlapped (point per decade)
            % Calculation of statistics
            %t=[0.0:dt:1800-dt];
            %x(:,1)=xx(:,i);
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
            psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % Check the Parseval theorem
%             display ('Test the Parseval Theorem')
%             display ('The value below should be near around 1.0 for the energy conservation')
%             sum(psd)./(n.*dt)
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
            %for js = 1:1:nn
    %         for js = 1:1:28
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
            % Plotting
            %loglog(f(:),fpsd(:,1),'xk',ax(:),fpsd_s(:,1),'-or')
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+50) = fpsd_s(ii,1);
            end
            clear ii fpds_s        

            % Spectra analysis
            x = temp(:,3);
            nx=4;   % the number of variables
            ii=3;    %
            dt=0.1;  % time difference between two conjecutive points (seconds)
            nst=16; % starting point of non-overlapping average
            pt_dec=10;  % how many points will be overlapped (point per decade)
            % Calculation of statistics
            %t=[0.0:dt:1800-dt];
            %x(:,1)=xx(:,i);
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
            psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % Check the Parseval theorem
%             display ('Test the Parseval Theorem')
%             display ('The value below should be near around 1.0 for the energy conservation')
%             sum(psd)./(n.*dt)
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
            %for js = 1:1:nn
    %         for js = 1:1:28
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
            % Plotting
            %loglog(f(:),fpsd(:,1),'xk',ax(:),fpsd_s(:,1),'-or')
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+96) = fpsd_s(ii,1);
            end
            clear ii fpds_s                

            % Spectra analysis
            x = temp(:,4);
            nx=4;   % the number of variables
            ii=3;    %
            dt=0.1;  % time difference between two conjecutive points (seconds)
            nst=16; % starting point of non-overlapping average
            pt_dec=10;  % how many points will be overlapped (point per decade)
            % Calculation of statistics
            %t=[0.0:dt:1800-dt];
            %x(:,1)=xx(:,i);
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
            psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % Check the Parseval theorem
%             display ('Test the Parseval Theorem')
%             display ('The value below should be near around 1.0 for the energy conservation')
%             sum(psd)./(n.*dt)
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
            %for js = 1:1:nn
    %         for js = 1:1:28
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
            % Plotting
            %loglog(f(:),fpsd(:,1),'xk',ax(:),fpsd_s(:,1),'-or')
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+142) = fpsd_s(ii,1);
            end
            clear ii fpds_s                

            % Spectra analysis
            x = temp(:,5);
            nx=4;   % the number of variables
            ii=3;    %
            dt=0.1;  % time difference between two conjecutive points (seconds)
            nst=16; % starting point of non-overlapping average
            pt_dec=10;  % how many points will be overlapped (point per decade)
            % Calculation of statistics
            %t=[0.0:dt:1800-dt];
            %x(:,1)=xx(:,i);
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
            psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
            f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
            for j=1:1:n2-1
                f(j) = real(j)/(real(n).*dt);
                psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            end
            f(n2) = real(n2)/(real(n).*dt);
            psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
            % Check the Parseval theorem
%             display ('Test the Parseval Theorem')
%             display ('The value below should be near around 1.0 for the energy conservation')
%             sum(psd)./(n.*dt)
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
            %for js = 1:1:nn
    %         for js = 1:1:28
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
            % Plotting
            %loglog(f(:),fpsd(:,1),'xk',ax(:),fpsd_s(:,1),'-or')
            clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi ii index j js k lf_ma
            clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
            for ii = 1:46
                result(i,ii+188) = fpsd_s(ii,1);
            end
            clear ii fpds_s                        

                fd2x = zeros(48,1);
    fd2y = zeros(48,1);

    for id2 = 1:1%48
    %     x = d300(((id2-1)*18000+1):id2*18000,3);
    %     y = d300(((id2-1)*18000+1):id2*18000,1);
        x = temp(:,3);
        y = temp(:,5);

        % Program History
        % NOV/30/2009   Original code (Jinkyu Hong)
        % Parameter setting
        nx=4;   % the number of variables
%         i=3;    %
        dt=0.1;  % time difference between two conjecutive points (seconds)
        nst=16; % starting point of non-overlapping average
        pt_dec=10;  % how many points will be overlapped (point per decade)
            % see Kaimal and Finnigan (1994) for the details
        % Calculation of statistics
        %t=[0.0:dt:1800-dt];
        %x(:,1)=xx(:,i);
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
        psd(:,1) = zeros(n2,1); % zero padding: power spectrum density

        copsd(:,1) = zeros(n2,1);

        f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
    %     h = waitbar(0,'1st. Please wait...'); % waitbar : je-woo
        for j=1:1:n2-1
    %         waitbar(j/n2) % waitbar : je-woo
            f(j) = real(j)/(real(n).*dt);
            psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);
            copsd(j,1) = (2.*n.*dt)*co(j+1,1)./coxy;
        end
    %     close(h) % waitbar : je-woo
        f(n2) = real(n2)/(real(n).*dt);
        psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
        copsd(n2,1) = (n.*dt)*co(n2+1,1)./coxy;

        % Check the Parseval theorem
%         display ('Test the Parseval Theorem')
%         display ('The value below should be near around 1.0 for the energy conservation')
%         sum(psd)./(n.*dt)
%         sum(copsd)./(n.*dt)    
        % frequency weighted power spectrum
        fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
        fpsd(:,1) = f(:).*psd(:,1);

        cofpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
        cofpsd(:,1) = f(:).*copsd(:,1);    
        % spectrum smoothing (high frequency band)
        logx=zeros(n2,1); % je-woo 2012-11-17
    %     h = waitbar(0,'2nd. Please wait...'); % waitbar : je-woo
        for j = nst:1:n2
    %         waitbar(j/n2) % waitbar : je-woo
            logx(j) = log10(f(j));
        end
    %     close(h) % waitbar : je-woo
        lf_mi = logx(nst);
        lf_ma = logx(n2);
        nn = int8((lf_ma - lf_mi).*real(pt_dec));
        ax=zeros(n2,1); % je-woo 2012-11-17
    %     h = waitbar(0,'3rd. Please wait...'); % waitbar : je-woo
        for js = 1:1:nst-1
    %         waitbar(js/nst) % waitbar : je-woo
            ax(js) = f(js);
        end
    %     close(h) % waitbar : je-woo
        fpsd_s=zeros(n2,1); % je-woo 2012-11-17
        cofpsd_s=zeros(n2,1); % je-woo 2012-11-17    
    %     h = waitbar(0,'4th. Please wait...'); % waitbar : je-woo
        for js = 1:1:nst-1
            fpsd_s(js,1) = fpsd(js,1);
            cofpsd_s(js,1) = cofpsd(js,1);
    %         waitbar(js/nst) % waitbar : je-woo
        end
    %     close(h) % waitbar : je-woo
        p = nst - 1;
        k = nst - 1;
        sp = 0.0;
    %     h = waitbar(0,'5th. Please wait...'); % waitbar : je-woo
        %for js = 1:1:nn
        %for js = 1:1:28
        for js = 1:1:45
    %         waitbar(js/45) % waitbar : je-woo
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
    %     close(h) % waitbar : je-woo
        clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi index j js k lf_ma
        clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar x
        clear y ym yp ystd yvar fyr fyi fy coxy cosp copsd cofpsd co

    %     hold on
    %     loglog(ax(1:46),fpsd_s(1:46,1),'-or')

        for id3 = 1:48
            fd2x(id3,id2) = ax(id3,1);
            fd2y(id3,id2) = fpsd_s(id3,1);
            cofd2y(id3,id2) = cofpsd_s(id3,1);
        end
        clear id3 ax fpsd_s cofpsd_s 
    end
    clear id2
    for jw = 1:46
        result(i,jw+234) = cofd2y(jw,1);
    end
    clear jw
% end
clear jewoo
            
        end
    end
%     close(wait);
%     clear i
%     clear wait
%     clear po_u po_v po_w po_Ts po_CO2
%     clear size_n size_var sonic_ang num_30min_n
    %clear temp temp_raw k n des_window des_sigma des_n_temp des_n_before

    %%
    % d2_11 = d2(828001:829800,:);
    % d60_11 = d60(828001:829800,:);
    % d140_11 = d140(828001:829800,:);
    % d300_11 = d300(828001:829800,:);
    % 
    % d2_12 = d2(900001:901800,:);
    % d60_12 = d60(900001:901800,:);
    % d140_12 = d140(900001:901800,:);
    % d300_12 = d300(900001:901800,:);
    % 
    % wtc(d2_11(:,8),d60_11(:,10))
    % wtc(d2_12(:,8),d60_12(:,10))
    % 
    % d2(:,5:11)=[];
    % d60(:,7:11)=[];
    % d140(:,5)=[];
    % d300(:,7:11)=[];




total_result{nlm,1}= dataName(nlm).name;
total_result{nlm,2}= result;
clear data result
end
clear nlm


