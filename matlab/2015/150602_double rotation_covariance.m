tic()
temp = data(350001:368000,:);
u_bar = mean(temp(:,5))
v_bar = mean(temp(:,6))
w_bar = mean(temp(:,7))

% double rotation
temp_wind = zeros(18000,3);
alpha = atan2(v_bar,u_bar);
beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.5)));
for j = 1:18000
    temp_wind(j,1) = cos(beta)*(cos(alpha)*temp(j,5)+sin(alpha)*temp(j,6))+sin(beta)*temp(j,7);
    temp_wind(j,2) = -sin(alpha)*temp(j,5)+cos(alpha)*temp(j,6);
    temp_wind(j,3) = -sin(beta)*(cos(alpha)*temp(j,5)+sin(alpha)*temp(j,6))+cos(beta)*temp(j,7);
%    temp(j,5) = temp_wind(j,5);
%    temp(j,6) = temp_wind(j,6);
%    temp(j,7) = temp_wind(j,7);
end
%clear u_bar v_bar w_bar alpha beta j temp_wind

u_bar2 = mean(temp_wind(:,1))
v_bar2 = mean(temp_wind(:,2))
w_bar2 = mean(temp_wind(:,3))


%% u'w'
uw=0;
for i = 1:18000
    uw = uw + (temp_wind(i,1)-u_bar2)*(temp_wind(i,3)-w_bar2);
end
uw=uw/18000
clear i
toc()
    

