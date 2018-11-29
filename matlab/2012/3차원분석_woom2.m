[n v] = size(data);
clear v

grid_n = 33;

x_min = min(data(:,1));
x_max = max(data(:,1));
x_ra = (x_max - x_min)/grid_n;
y_min = min(data(:,2));
y_max = max(data(:,2));
y_ra = (y_max - y_min)/grid_n;

result = zeros(grid_n, grid_n);
for i = 1:n
    x = (data(i,1)-x_min)/x_ra - mod((data(i,1)-x_min)/x_ra,1);
    y = (data(i,2)-y_min)/y_ra - mod((data(i,2)-y_min)/y_ra,1);
    if x<=0
        x=1;
    end
    if y<=0
        y=1;
    end
    result(y,x) = result(y,x)+1;
end
for i = 1:grid_n
    for j = 1:grid_n
        result(i,j) = result(i,j)*100/n;
    end
end
contourf(result)



