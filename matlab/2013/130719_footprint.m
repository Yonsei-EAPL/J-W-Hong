size_x = 1529;

%% constant 
g = 9.799364;
k = 0.4;
z0 = 0.6;
zd = 4;
zm = 10;
z= zm-zd;
zu = z*(log(z/z0)-1.0+(z0/z));

%% extract
temp = zeros(4001,4001);

% i = 1;
for i = 1:1529
    mean_u = data(i,1);
    std_v = data(i,2);
    wd = data(i,3);
    wd = wd/180*pi();
    u_star = data(i,4);
    H = data(i,5);
    if data(i,6)<100
        mean_T = data(i,6)+273.15;
    else
        mean_T = data(i,6);
    end

    L = -1*u_star^3*mean_T/k/g/(H/1230);
    zoL = zu/L;
    if zoL<-0.04
        D=0.28;
        P=0.59;
    elseif zoL>0.04
        D=2.44;
        P=1.33;
    else
        D=0.97;
        P=1.0;
    end

    cons = 1/k^2*D*zu^P*abs(L)^(1-P);
    x_max = cons/2;

    for x = 1:2000
        fx = 1/x^2*cons*exp(-1/x*cons);
        a=mean_u/(std_v*x);
%         thre=abs(3/a);
        for y=1:2000
%             if y<thre
                fy=a/((2*pi())^(0.5))*exp(-1/2*(y*a)^2);
                f=fx*fy;
                if f>0
                    x1=x;
                    y1=y;
                    xx=y1*cos(wd)-x1*sin(wd);
                    yy=y1*sin(wd)-x1*cos(wd);
                    xx2=-1*y1*cos(wd)-x1*sin(wd);
                    yy2=-1*y1*sin(wd)-x1*cos(wd);
                    xx=xx-mod(xx,1)+2000;
                    yy=yy-mod(yy,1)+2000;
                    xx2=xx2-mod(xx2,1)+2000;
                    yy2=yy2-mod(yy2,1)+2000;
                    if ((xx<=4001)&&(xx>0))&&((yy<=4001)&&(yy>0))
                        temp(xx,yy) = temp(xx,yy)+f;
                    end
                    if ((xx2<4001)&&(xx2>0))&&((yy2<4001)&&(yy2>0))
                        temp(xx2,yy2) = temp(xx2,yy2)+f;
                    end
                end
%             end
        end
    end
    
end