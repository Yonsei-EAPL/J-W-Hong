%% history
% 2013-05-22 J-W Hong, for Eunpyeong Newtown

%% double rotation 
% for 10Hz, 30min
% data(18000,3) ; 1-u 2-v, 3-w

u_bar = mean(data(:,1));
v_bar = mean(data(:,2));
w_bar = mean(data(:,3));

alpha = atan(v_bar/u_bar);
beta = atan(w_bar/(sqrt(u_bar^2+v_bar^2)));

rotated_uvw=zeros(18000,3);

for i = 1:18000
	rotated_uvw(i,1) = cos(beta)*(cos(alpha)*data(i,1)+sin(alpha)*data(i,2))+sin(beta)*data(i,3);
	rotated_uvw(i,2) = -sin(alpha)*data(i,1)+cos(alpha)*data(i,2);
	rotated_uvw(i,3) = -sin(beta)*(cos(alpha)*data(i,1)+sin(alpha)*data(i,2))+cos(beta)*data(i,3);
end


%% Wind-direction 

wd_sonic = 220; % sonic direction
wd = zeros(18000,1);
for i = 1:18000
   if data(i,1) >0
       if data(i,2) >0 % 4th quadrant
           wd(i,1)=wd_sonic-1*atan(data(i,2)/data(i,1))/pi()*180;
       elseif data(i,2) <0 % 1st quadrant
           wd(i,1)=wd_sonic+atan((-1*data(i,2))/data(i,1))/pi()*180;
       end
   elseif data(i,1) <0
       if data(i,2) >0 % 3rd quadrant
           wd(i,1)=wd_sonic-180+atan(data(i,2)/(-1*data(i,1)))/pi()*180;
       elseif data(i,2) <0 % 2nd quadrant 
           wd(i,1)=wd_sonic-180-atan((-1*data(i,2))/(-1*data(i,1)))/pi()*180;
       end       
   end
end
for i = 1:18000
    if wd(i,1) >360
        wd(i,1) = wd(i,1) -360;
    elseif wd(i,1) <0
        wd(i,1) = wd(i,1) +360;
    end
end
clear i