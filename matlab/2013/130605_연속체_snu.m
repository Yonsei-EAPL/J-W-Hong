%%

L=1;
alpha = 3.5*10^(-7);
T0 = 283.15;
dx = 0.01;
t = 1000000;
dt = 100;

rga = 0;
rgb = L/dx;
D = alpha*dt/dx^2;
m = t/dt;

temp = zeros(rgb+1,2);

%%

temp(1,1) = T0;

for i = 2:rgb+1
    temp(i,1) = 0;
    temp(i,2) = 0;
end

for i = 1:m
    for j = 2:rgb
        temp(j,2) = D*temp(j+1,1)+(1-2*D)*temp(j,1)+D*temp(j-1,1);
    end
    temp(rgb,2) = temp(rgb-1,2);
    temp(1,2) = T0;
    for j = rga+1:rgb
        temp(j,1) = temp(j,2);
    end
    plot(temp(:,1))
    hold on
end

