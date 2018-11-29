n_size = 270;

paramEsts = gevfit(x);
ka = paramEsts(1,1)
si = paramEsts(1,2)
mu = paramEsts(1,3)
clear paramEsts

Fx = zeros(n_size,1);
for i = 1:n_size
    Fx(i,1) = exp(-(1+ka*(x(i,1)-mu)/si)^(-ka^(-1)));
end
clear i
