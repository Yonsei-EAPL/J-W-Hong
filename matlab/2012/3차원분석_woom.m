temp = [2574];
grid_n = 33;

[n v] = size(data);
n_v = v/2;
if n_v>5
    n_v = 5;
end
for i = 1:n_v
    open('h1(4).fig')
end
clear n v i

y_min = -8;
y_max = 8;
y_ra = (y_max - y_min)/grid_n;
x_min = 40;
x_max = 80;
x_ra = (x_max - x_min)/grid_n;

for i = 1:n_v
    result = zeros(grid_n,grid_n);
    for j = 1:temp(i,1)
        x = (data(j,(i-1)*2+1)-x_min)/x_ra - mod((data(j,(i-1)*2+1)-x_min)/x_ra,1);
        y = (data(j,(i-1)*2+2)-y_min)/y_ra - mod((data(j,(i-1)*2+2)-y_min)/y_ra,1);
        if x <=0
            x = 1;
        elseif x>33
            x = 33;
        end
        if y <=0
            y = 1;
        elseif y>33
            y=33;
        end
        result(y,x) = result(y,x)+1;
    end
    for j = 1:grid_n
        for k = 1:grid_n
            result(j,k) = result(j,k)*100/temp(i,1);
        end
    end
    max(max(result))
    figure(i)
    pause()
    hold on
    contourf(result)
end
clear i j k x y n y_min y_max x_min x_max y_ra x_ra result


