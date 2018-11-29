dataDir = 'E:\EAPL\ep\_csv'; % folder name % EP
dataName = dir(fullfile(dataDir, 'CSV_3302.ts_data_2013_*'))'; % EP 1st system

for t=1:length(dataName) %n_st:n_en% 여러개의 CSV 파일을 연속적으로 처리하기 위한 loop
    
    data = importdata(fullfile(dataDir, dataName(t).name));  
    co2 = data(:,10);
    if length(co2)==864000
        for iii=1:48
            x=co2((iii-1)*18000+1:iii*18000,:);

        % !!! before run this code, move your data to varable x !!! (ex. x=data(:,1);)
        % dt means time-interval of your raw-data (ex. if 20Hz data then, dt=0.05;)

        %% Program History
        % NOV/30/2009   Original code (Jinkyu Hong)
        %clear all
        %clc
        %clf

        %%
        % Reading input file
        %xx=dlmread('input.dat');


        %%
        % Parameter setting
        nx=4;   % the number of variables
        i=3;    %
        dt=0.1;  % time difference between two conjecutive points (seconds)
        %dt=60;

        nst=16; % starting point of non-overlapping average
        pt_dec=10;  % how many points will be overlapped (point per decade)
            % see Kaimal and Finnigan (1994) for the details


        %%
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

        %%
        % FFT
        fx(:,1)=fft(xp(:,1));
        fxr(:,1)=real(fx(:,1))./real(n);
        fxi(:,1)=imag(fx(:,1))./real(n);

        %%
        % Power spectrum 
        s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;


        %%
        % Power spectrum density
        psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
        f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)

        % h = waitbar(0,'1st. Please wait...'); % waitbar : je-woo
        for j=1:1:n2-1
        %     waitbar(j/n2) % waitbar : je-woo
            f(j) = real(j)/(real(n).*dt);
            psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
        end
        % close(h) % waitbar : je-woo

        f(n2) = real(n2)/(real(n).*dt);
        psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);

        % Check the Parseval theorem
        display ('Test the Parseval Theorem')
        display ('The value below should be near around 1.0 for the energy conservation')
        sum(psd)./(n.*dt)

        % frequency weighted power spectrum
        fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
        fpsd(:,1) = f(:).*psd(:,1);

        %%
        % spectrum smoothing (high frequency band)

        logx=zeros(n2,1); % je-woo 2012-11-17
        % h = waitbar(0,'2nd. Please wait...'); % waitbar : je-woo
        for j = nst:1:n2
        %     waitbar(j/n2) % waitbar : je-woo
            logx(j) = log10(f(j));
        end
        % close(h) % waitbar : je-woo

        lf_mi = logx(nst);
        lf_ma = logx(n2);
        nn = int8((lf_ma - lf_mi).*real(pt_dec));

        ax=zeros(n2,1); % je-woo 2012-11-17
        % h = waitbar(0,'3rd. Please wait...'); % waitbar : je-woo
        for js = 1:1:nst-1
        %     waitbar(js/nst) % waitbar : je-woo
            ax(js) = f(js);
        end
        % close(h) % waitbar : je-woo

        fpsd_s=zeros(n2,1); % je-woo 2012-11-17
        % h = waitbar(0,'4th. Please wait...'); % waitbar : je-woo
        for js = 1:1:nst-1
            fpsd_s(js,1) = fpsd(js,1);
        %     waitbar(js/nst) % waitbar : je-woo
        end
        % close(h) % waitbar : je-woo

        p = nst - 1;
        k = nst - 1;
        sp = 0.0;


        % h = waitbar(0,'5th. Please wait...'); % waitbar : je-woo
        %for js = 1:1:nn
        for js = 1:1:28
        %for js = 1:1:45

        %     waitbar(js/45) % waitbar : je-woo

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
        % close(h) % waitbar : je-woo
        %%
        % Plotting
        %figure()
        %loglog(f(:),fpsd(:,1),'xk',ax(:),fpsd_s(:,1),'-or')
        if iii==1
            figure()
        end
        loglog(ax(:),fpsd_s(:,1),'-o','color',rand(1,3))
        hold on

        clear ans dt flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
        clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x

        if iii==48
            y=[0.0100000000000000,1;10,0.0100000000000000;];
            hold on
            plot(y(:,1),y(:,2))
        end

        %loglog(f(:),fpsd2pi(:,1),'xk',ax(:),fpsd_s2pi(:,1),'-or')

        end
    end
    title(sprintf('%s',dataName(t).name));    
    saveas(gcf,sprintf('%s.png',dataName(t).name));    

    
    
    close
end
clear all
clc

