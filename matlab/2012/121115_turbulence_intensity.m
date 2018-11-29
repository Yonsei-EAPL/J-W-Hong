
%% program for "turbulence intensity" analysis
%   2012-11-15, Je-Woo Hong (advised, Prof.Jinkyu Hong)

%% information : change yourself

range = 24; % data range (unit:hr)
period = 5; % averaging period for one point (unit:min)
observ = 20; % observ freq. (unit:Hz)

% Also, Check the position of your rawdata(wind and sonic-temperature), 
% and Change the "pick data-set" in "Calculation" part.

%   u : along-wind speed,
%   v : cross-wind speed, 
%   w : vertical-wind speed,
%   ts (or t) : sonic temperature,

%% variables

no_point = range * 60 / period; % number of points for analysis range
no_data = period * 60 * observ; % number of data for one point
range = range * 60 * 60 * observ;
period = period * 60;

mean_u = zeros(no_point,1);
mean_t = zeros(no_point,1);
mean_w = 0;
sigma_w = zeros(no_point,1); % sigma means standard deviation
sigma_t = zeros(no_point,1);

temp_u = zeros(no_data,1);
temp_v = zeros(no_data,1);
temp_w = zeros(no_data,1);
temp_t = zeros(no_data,1);
temp = zeros(no_data,1);

%% calclation

for i = 1:no_point
    
    % pick data-set
    for k = 1:no_data
        temp_u(k,1) = raw_data((i-1)*no_data+k,1); % !! change, if your raw data is different
        temp_v(k,1) = raw_data((i-1)*no_data+k,2); % !! change, if your raw data is different
        temp_w(k,1) = raw_data((i-1)*no_data+k,3); % !! change, if your raw data is different
        temp_t(k,1) = raw_data((i-1)*no_data+k,4); % !! change, if your raw data is different
    end
   
    % double rotation
    u_bar = mean(temp_u);
    v_bar = mean(temp_v);
    w_bar = mean(temp_w);
    alpha = atan(v_bar/u_bar);
    beta = atan(w_bar/(sqrt(u_bar^2+v_bar^2)));
    for k = 1:no_data
        temp_u(k,1) = cos(beta)*(cos(alpha)*temp_u(k,1)+sin(alpha)*temp_v(k,1))+sin(beta)*temp_w(k,1);
        temp_v(k,1) = -sin(alpha)*temp_u(k,1)+cos(alpha)*temp_v(k,1);
        temp_w(k,1) = -sin(beta)*(cos(alpha)*temp_u(k,1)+sin(alpha)*temp_v(k,1))+cos(beta)*temp_w(k,1);
    end
    
    % mean wind velocity
    for k = 1:no_data
        temp(k,1) = sqrt(temp_u(k,1)^2+temp_v(k,1)^2+temp_w(k,1)^2);
    end
    mean_u(i,1) = mean(temp);
    temp = zeros(no_data,1);
    
    % sigma(standard deviation) w
    mean_w = mean(temp_w);
    for k = 1:no_data
        temp(k,1) = (temp_w(k,1) - mean_w)^2;
    end
    sigma_w(i,1) = sqrt(mean(temp));
    temp = zeros(no_data,1);
    mean_w = 0;
    
    % mean temperature
    mean_t(i,1) = mean(temp_t);
    
    % sigma t
    for k = 1:no_data
        temp(k,1) = (temp_t(k,1) - mean_t(i,1))^2;
    end
    sigma_t(i,1) = sqrt(mean(temp));
    temp = zeros(no_data,1);
        
end

% turbulence intensity
i_w = zeros(no_point,2);
i_t = zeros(no_point,2);
for i = 1:no_point
    i_w(i,1)=(i-1)*period+period*0.5; % time for x-axis
    i_w(i,2)=abs(sigma_w(i,1)/mean_u(i,1)); % turbulence intensity of w 
    i_t(i,1)=(i-1)*period+period*0.5; % time for x-axis
    i_t(i,2)=abs(sigma_t(i,1)/mean_t(i,1)); % turbulence intensity of ts 
end

%% plot
subplot(2,1,1)
plot(i_w(:,1),i_w(:,2)); % subplot(up) : w
subplot(2,1,2)
plot(i_t(:,1),i_t(:,2)); % subplot(down) : ts

%% clear
clear i k u_bar v_bar w_bar alpha beta temp temp_u temp_v temp_w temp_t mean_w sigma_w sigma_t mean_u mean_t no_point no_data period range
	