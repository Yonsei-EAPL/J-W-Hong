%% footprint _ EAPL
% J-W Hong and prof. J.Hong modified from Keunmin Lee ans prof. J.Hong 's fortran code
%   based on, Hsieh et al.(2000) and Laubach ans Kelliher(2004)
% 2014-08-14, estimate land_fraction (using LULC figure)
    % 1 : building
    % 2 : roads
    % 3 : other impervious
    % 4 : water
    % 5 : forest (tree)
    % 6 : grass
    % 7 : garden
    % 8 : bare soil


%% information
length = 659;


%% constant
g=9.799364; % gravity acceleration from KIST
k=0.4; % von Karman constant


%% bin for result
fxy_final = zeros(4001,4001); % final result rotated ans accumulated
LULC_result = zeros(length,9); % 1to 8 and out of boundary


%% calculation
wait = waitbar(0,'...Footprint Analysis...');

for i = 1:length
    waitbar((i/length),wait,sprintf('%3f',i))

    % before rotation
    fxy = zeros(2001,2001); % result from each time-step
    
    u = abs(data(i,1));
    std_v = data(i,2);
    wd = data(i,3);
    u_star = data(i,4);
    H = data(i,5);
    T = data(i,6); % degree C

%     z0 = 0.7; % SF
%     zd = 4.9; % SF
%     zm = 10; % SF
    z0 = 0.5; % EP
    zd = 5.2; % EP
    zm = 30; % EP
    z = zm-zd;
    zu = z*(log(z/z0)-1+(z0/z));
    
    L = -u_star^3*(T+273.15)/k/g/H*1230;
    zoL = zu/L;
    if zoL>=0.04
        D = 2.44;
        P = 1.33;
    elseif zoL<=-0.04
        D = 0.28;
        P = 0.59;
    else
        D = 0.97;
        P = 1.0;
    end
    
    cons = 1/k^2*D*zu^P*abs(L)^(1-P);
    mx = 1/2*cons;
    
    fx = zeros(2000,1);
    for x = 1:2000
        fx(x,1) = 1/x^2*cons*exp(-1/x*cons);
        a=u/std_v/x;
        thre = 4/a;
        if thre<0
            thre = thre-mod(thre,1);
            if abs(thre)>1000
                thre = -1000;
            end
            if thre== 0
                thre = -1;
            end
            for y=thre:-thre
                fy = a/(2*pi())^(0.5)*exp(-1/2*(abs(y)*a)^2);
                f= fx(x,1)*fy;
                fxy(x,1001+y) = f;
            end
        elseif thre>0
            thre = thre-mod(thre,1);
            if abs(thre)>1000
                thre = 1000;
            end            
            if thre == 0
                thre = 1;
            end
            for y=-thre:thre
                fy = a/(2*pi())^(0.5)*exp(-1/2*(abs(y)*a)^2);
                f= fx(x,1)*fy;
                fxy(x,1001+y) = f;
            end
        end
    end
    
    % rotation ans accumulation
    wd = wd/180*pi();    
    for x = 1:2001
        for y = 1:2001
            if fxy(x,y)==0
            else
                x1 = x;
                y1 = y-1001;
                x2 = cos(wd)*y1 + sin(wd)*x1;
                y2 = -sin(wd)*y1 + cos(wd)*x1;
                x2 = x2 - mod(x2,1);
                if x2<=-2000
                    x2 = -2000;
                elseif x2>=2000
                    x2 = 2000;
                end
                y2 = y2 - mod(y2,1);
                if y2<=-2000
                    y2 = -2000;
                elseif y2>=2000
                    y2 = 2000;
                end
                fxy_final(2001+x2,2001+y2) = fxy_final(2001+x2,2001+y2)+fxy(x,y);
                
                if (abs(x2)<=499)&&(abs(y2)<=499)
                    LULC_result(i,LULC(500+x2,500+y2)) = LULC_result(i,LULC(500+x2,500+y2))+fxy(x,y);
                end
            end
        end
    end
    temp_sum = 0;
    for j = 1:8
        temp_sum = temp_sum + LULC_result(i,j);
    end
    LULC_result(i,9) = 1 - temp_sum;
%     figure(1)
%     hold on
%     plot(fx);
end
clear i j k temp_sum D H L P T a cons f fx fxy fy g mx std_v u u_star wain wd x x1 x2 y y1 y2 z z0 zd zm zoL zu
temp_fxy = fxy_final;
fxy_final = zeros(4001,4001);
for i = 1:4001
    for j = 1:4001
        if temp_fxy(j,i)>0
            fxy_final(i,j) = temp_fxy(j,i)/length*100;
        end
    end
end
clear i j temp_fxy
close(wait);
% temp = 0;
% for i = 1:4001
%     for j = 1:4001
%         if fxy_final(i,j)>=0.001
%             temp = temp + fxy_final(i,j);
%         end
%     end
% end
% clear i j
