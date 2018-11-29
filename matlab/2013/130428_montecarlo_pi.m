
l=10000;
result = zeros(l,1);

for j = 1:l

    inside = 0;
    outside = 0;

    for i = 1:l

        x= rand();
        y= rand();

        if x^2+y^2 <= 1
            inside = inside+1;
        else 
            outside = outside + 1;
        end

    end
    result(j,1) = inside/l*4;

end

mean(result)
