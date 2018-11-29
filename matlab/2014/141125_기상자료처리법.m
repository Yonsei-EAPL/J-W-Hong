%% History
% 2014-11-25, Je-Woo Hong, Spectrum Analysis 

%% data : y
dt = 1800; % data interval in second
y_bar = mean(y); % y0
[n, m] = size(y); % total data length information
clear m
if mod(n,2)==0
    k = n/2; % even
else
    k = (n-1)/2; % odd
end
a_k = zeros(k,1); 
b_k = zeros(k,1);
c_k = zeros(k,1);
phi_k = zeros(k,1); % phase
r_k_2 = zeros(k,1); % R square
s_k_2 = zeros(k,1);% variance

%% main run
for i = 1:k
    for j = 1:n
        a_k(i,1) = a_k(i,1)+ y(j,1)*cos(2*pi()*i*j/n);
        b_k(i,1) = b_k(i,1)+ y(j,1)*sin(2*pi()*i*j/n);
    end
    a_k(i,1) = a_k(i,1)*2/n;
    b_k(i,1) = b_k(i,1)*2/n;
    c_k(i,1) = (a_k(i,1)^2 + b_k(i,1)^2)^(0.5);
    r_k_2(i,1) = (a_k(i,1)^2 + b_k(i,1)^2);
    if a_k(i,1)>0
        phi_k(i,1) = atan(b_k(i,1)/a_k(i,1));
    elseif a_k(i,1)<0
        phi_k(i,1) = atan(b_k(i,1)/a_k(i,1))+ pi();
    else
        pi()/2;
    end
    if phi_k(i,1)>2*pi()
        phi_k(i,1) = phi_k(i,1) - 2*pi();
    elseif phi_k(i,1)<0
        phi_k(i,1) = phi_k(i,1) + 2*pi();
    end
    if i ==k
        s_k_2(i,1) = a_k(i,1)^2;
    else
        s_k_2(i,1) = (a_k(i,1)^2 + b_k(i,1)^2)/2;
    end
end
clear i j 
total_variance = sum(s_k_2);

%% frequency 
f_k = zeros(k,1);
for i = 1:k
    f_k(i,1) = (1/((n/i)*dt));
end
clear i

%% plotting (log scale)
figure()
loglog(f_k(:,1),c_k(:,1),'-xk')
