xi = surface_data(864001:936000,5:7);
x = zeros(72000,1);
for i = 1:72000
    x(i,1) = (xi(i,1)^2 + xi(i,2)^2 + xi(i,3)^2)^(0.5);
end
clear i xi

xi = m60_data(864001:936000,5:7);
x = zeros(72000,1);
for i = 1:72000
    x(i,1) = (xi(i,1)^2 + xi(i,2)^2 + xi(i,3)^2)^(0.5);
end
clear i xi

xi = m140_data(864001:936000,5:7);
x = zeros(72000,1);
for i = 1:72000
    x(i,1) = (xi(i,1)^2 + xi(i,2)^2 + xi(i,3)^2)^(0.5);
end
clear i xi

xi = m300_data(864001:936000,5:7);
x = zeros(72000,1);
for i = 1:72000
    x(i,1) = (xi(i,1)^2 + xi(i,2)^2 + xi(i,3)^2)^(0.5);
end
clear i xi

xi = surface_data(864001:936000,5:7);
x = zeros(72000,1);
for i = 1:72000
    x(i,1) = (xi(i,1)^2 + xi(i,2)^2 + xi(i,3)^2)^(0.5);
end
clear i xi

xi = m60_data(864001:936000,5:7);
x = zeros(72000,1);
for i = 1:72000
    x(i,1) = (xi(i,1)^2 + xi(i,2)^2 + xi(i,3)^2)^(0.5);
end
clear i xi

xi = m140_data(864001:936000,5:7);
x = zeros(72000,1);
for i = 1:72000
    x(i,1) = (xi(i,1)^2 + xi(i,2)^2 + xi(i,3)^2)^(0.5);
end
clear i xi

xi = m300_data(864001:936000,5:7);
x = zeros(72000,1);
for i = 1:72000
    x(i,1) = (xi(i,1)^2 + xi(i,2)^2 + xi(i,3)^2)^(0.5);
end
clear i xi






open('spectrum comparision.fig')
hold on
loglog(surface_ax(:),surface_fpsd_s(:,1),'-or')
hold on
loglog(m60_ax(:),m60_fpsd_s(:,1),'-om')
hold on
loglog(m140_ax(:),m140_fpsd_s(:,1),'-og')
hold on
loglog(m300_ax(:),m300_fpsd_s(:,1),'-ob')