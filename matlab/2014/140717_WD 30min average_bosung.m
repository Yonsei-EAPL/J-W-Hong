%% History
% 140717 Mr.Keunmin Lee; code

%% Average Wind direction
% input : WS, WD
% u,v : u-vector, v-vector
u=zeros(length(WD),2);
v=zeros(length(WD),2);
[size_x, size_y]=size(WS);

for j = 1:size_y
    for i = 1:size_x
        u(i,j)=WS(i,j)*-1*sin(WD(i)*pi/180);
        v(i,j)=WS(i,j)*-1*cos(WD(i)*pi/180);
    end
    clear i

    %30min average
    x=0;
    for i = 1:size_x
        if mod(i,30)==1
            x=x+1;
            y=1;
            u_temp(x,y)= u(i,j);
            v_temp(x,y)= v(i,j);
        else
            y=y+1;
            u_temp(x,y)= u(i,j);
            v_temp(x,y)= v(i,j);
        end
    end
    clear i x y

    u_bar(:,j)=mean(u_temp,2);
    v_bar(:,j)=mean(v_temp,2);
    clear u_temp v_temp
end

clear j size_x size_y u v




