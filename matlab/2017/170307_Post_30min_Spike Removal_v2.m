%%  Input
% input : data(time,Fc,Qh,Qe,Kdown)
temp(:,1) = data(:,2); % Fc
temp(:,2) = data(:,3); % Qh
temp(:,3) = data(:,4); % Qe
kdown(:,1) = data(:,5); % Kdown

%[x,y] = find(temp==0); temp(x,y)=NaN;
[x] = find(temp==0); temp(x)=NaN;
clear x 
[x] = find(kdown==0); kdown(x)=NaN;
clear x y 
[nd nv] = size(temp);



%% Define Constants and Physical Ranges for Spike Detection

thrsh_Fc = 4.0;         % threshold for removing Fc spike using median
thrsh_Qh = 4.0;         % threshold for removing Qh spike using median
thrsh_Qe = 4.0;         % threshold for removing Qe spike using median

nday = 7;  %21 original            % processing period as the # of days
nod = 48;               % data number per day
d_range = nday*nod;     % processing range as the # of time
clear nday nod

kdown_limit = 5;         % downward radiation for night,daytime separation (W/m2)

%  Physical Ranges --------------------------------------------------------

f_Fc_upp = 100; % upper limit of Fc (umol/m2/s) 
f_Fc_low = -50; % lower limit of Fc (umol/m2/s) 

f_Qh_upp = 500;  % upper limit of Qh (W/m2) 
f_Qh_low = -200;  % lower limit of Qh (W/m2) 

f_Qe_upp = 500;  % upper limit of Qe (W/m2) 
f_Qe_low = -150;  % lower limit of Qe (W/m2) 

%%

% Calculation diferences for spike removal
d = zeros(nd,nv);
for j = 1:nv
	for i = 2:nd-1
        d(i,j) = (temp(i,j)-temp(i-1,j)) - (temp(i+1,j)-temp(i,j));
    end
	d(1,j) = d(2,j); d(nd,j)=d(nd-1,j); % 양 끝값 처리
end
        
clear i j    
    
% Time index calculation for spike removal
time = data(:,1);
excel = datenum(1900,1,1,0,0,0);
mat_time(1:length(time),1) = time(1:length(time),1) + excel - 2; 
time_class = datevec(mat_time);

clear excel mat_time time 
    
time_index = zeros(nd,1); % time index (0 : day, 1 : night)
for i = 1:nd
	if kdown(i)==NaN;
        if (time_class(i,4)<7||time_class(i,4)>18)
            time_index(i) = 1;
        end
    elseif (kdown(i)<=kdown_limit)&&((time_class(i,4)<10)||(time_class(i,4)>16))
         time_index(i) = 1;
	end
end
clear i time_class

% 각 30분 데이터를 n_range 설정에 따라 b개로 나누어 Spike Detection하기 위한 작업
num_bin = zeros(ceil(nd/d_range),3); 
for x = 1:nd
	if (x==1)
        b = 1;
        num_bin(b,1) = b;
        num_bin(b,2) = 1;
        num_bin(b,3) = 1;
    else
        if (num_bin(b,2)==d_range)
            b = b+1;
            num_bin(b,1) = b;
            num_bin(b,2) = 1;
            num_bin(b,3) = num_bin(b-1,3)+1;
        else
            num_bin(b,2) = num_bin(b,2)+1;
            num_bin(b,3) = num_bin(b,3)+1;
        end
    end 
end
        
clear x

%------------------------------------------------------------
%                    Main Processing
%------------------------------------------------------------
% Physical Range Detection Processing
flag = zeros(nd,nv);    % Spike flag
NaN_num = zeros(1,nv);  % Spike number

for x = 1:b
	b2 = num_bin(x,2); % 각 loop당 데이터 개수
    b3 = num_bin(x,3); % 누적 데이터 개수
    b4 = b3-b2+1;      % 각 loop당 데이터 시작 시점
    
    t_total = time_index(b4:b3);  
    t_day = numel(find(t_total==0));   % 각 loop당 낮 데이터 개수
    t_night = numel(find(t_total==1)); % 각 loop당 밤 데이터 개수
    clear t_total
    
    % Physical range-------------------------------------
    for i = b4:b3
        % Fc
        if(temp(i,1)>=f_Fc_upp);  flag(i,1) = 1; end
        if(temp(i,1)<=f_Fc_low);  flag(i,1) = 1; end
        % Qh
        if(temp(i,2)>=f_Qh_upp);  flag(i,2) = 1; end 
        if(temp(i,2)<=f_Qh_low);  flag(i,2) = 1; end                     
        % Qe
        if(temp(i,3)>=f_Qe_upp);  flag(i,3) = 1; end                         
        if(temp(i,3)<=f_Qe_low);  flag(i,3) = 1; end
    end
        
    % Time index ---------------------------------------
    for k = 1:2 
    if k==1 % Daytime
        dd = zeros(t_day,nv); ddd = zeros(t_day,nv);
        n=1;
        for i = b4:b3
            if time_index(i)==0
                dd(n,:) = d(i,:);
                n = n+1;
            end
        end
        clear i n
        
        Md = nanmedian(dd);
        for y = 1:t_day
            ddd(y,1:nv) = dd(y,1:nv)-Md(1,:);
        end
        clear y
        MAD = nanmedian(abs(ddd));
        
        cr_Fc = thrsh_Fc/0.6745*MAD(1,1);
        cr_Qh = thrsh_Qh/0.6745*MAD(1,2);
        cr_Qe = thrsh_Qe/0.6745*MAD(1,3);        
            
        for j = 1:nv
            if j==1
                for i = b4:b3
                    if (time_index(i)==0)
                        if (d(i,j)<Md(1,j)-cr_Fc)||(d(i,j)>Md(1,j)+cr_Fc)
                            flag(i,j) = 1;
                        end       
                    end
                end
                    
            elseif j==2
                for i = b4:b3
                    if (time_index(i)==0)
                        if (d(i,j)<Md(1,j)-cr_Qh)||(d(i,j)>Md(1,j)+cr_Qh)
                            flag(i,j) = 1;
                        end       
                    end
                end
                
            else
                for i = b4:b3
                    if (time_index(i)==0)
                        if (d(i,j)<Md(1,j)-cr_Qe)||(d(i,j)>Md(1,j)+cr_Qe)
                            flag(i,j) = 1;
                        end       
                    end
                end
            end
        end
        
    	clear dd ddd Md MAD cr_Fc cr_Qh cr_Qe i j
        
        
    else  % Nighttime
        dd = zeros(t_night,nv); ddd = zeros(t_night,nv);
        n = 1;
        for i = b4:b3
            if time_index(i)==1
                dd(n,:) = d(i,:);
                n = n+1;
            end
        end
        clear i n
        
        Md = nanmedian(dd);
        for y = 1:t_night
            ddd(y,1:nv) = dd(y,1:nv)-Md(1,:);
        end
        clear y
        MAD = nanmedian(abs(ddd));
        
        cr_Fc = thrsh_Fc/0.6745*MAD(1,1);
        cr_Qh = thrsh_Qh/0.6745*MAD(1,2);
        cr_Qe = thrsh_Qe/0.6745*MAD(1,3);        

        for j = 1:nv
            if j==1
                for i = b4:b3
                    if (time_index(i)==1)
                        if (d(i,j)<Md(1,j)-cr_Fc)||(d(i,j)>Md(1,j)+cr_Fc)
                            flag(i,j) = 1;
                        end       
                    end
                end
                    
            elseif j==2
                for i = b4:b3
                    if (time_index(i)==1)
                        if (d(i,j)<Md(1,j)-cr_Qh)||(d(i,j)>Md(1,j)+cr_Qh)
                            flag(i,j) = 1;
                        end       
                    end
                end
                
            else
                for i = b4:b3
                    if (time_index(i)==1)
                        if (d(i,j)<Md(1,j)-cr_Qe)||(d(i,j)>Md(1,j)+cr_Qe)
                            flag(i,j) = 1;
                        end       
                    end
                end
            end
        end
    	clear dd ddd Md MAD cr_Fc cr_Qh cr_Qe i j
    end
    end
    clear k
    clear b2 b3 b4 t_day t_night       
    end    
    clear x b time_index
        
%-------------------------------------------------------------------------    
    % Spike Removal Processing
	for  j = 1:nv
        for i = 1:nd
            if (flag(i,j)==1)
                temp(i,j) = NaN;
                NaN_num(1,j) = NaN_num(1,j)+1;
            end
        end
    end
      
    clear i j 
    clear nv nd
    
    
    
    
%%

    %--------------------------------------------------------------------------
    figure(1)
    xlim = [1,18000];
    ylim = [-100,200];
    plot(data(:,2),'-xk')
    hold on
    plot(temp(:,1),'og')
    set(gca,'XLim',xlim(:))
    set(gca,'YLim',ylim(:))
    ylabel(' Fc(umol/m2/s)')
    xlabel('data point')
    echo off
    
        %--------------------------------------------------------------------------
    figure(2)
    xlim = [1,18000];
    ylim = [-200,500];
    plot(data(:,3),'-xk')
    hold on
    plot(temp(:,2),'or')
    set(gca,'XLim',xlim(:))
    set(gca,'YLim',ylim(:))
    ylabel(' Qh(W/m2)')
    xlabel('data point')
    echo off
    
        %--------------------------------------------------------------------------
    figure(3)
    xlim = [1,18000];
    ylim = [-150,500];
    plot(data(:,4),'-xk')
    hold on
    plot(temp(:,3),'ob')
    set(gca,'XLim',xlim(:))
    set(gca,'YLim',ylim(:))
    ylabel(' Qe(W/m2)')
    xlabel('data point')
    echo off
    