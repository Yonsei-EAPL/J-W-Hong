result_seoul = zeros(35,9);

[x y] = size(result2006);
for j = 1:35
    for i = 1:x
        if result2006(i,1) == aws_seoul(j,1)
            result_seoul(j,2) = result2006(i,2);
        end
    end
end

[x y] = size(result2007);
for j = 1:35
    for i = 1:x
        if result2007(i,1) == aws_seoul(j,1)
            result_seoul(j,4) = result2007(i,2);
        end
    end
end

[x y] = size(result2008);
for j = 1:35
    for i = 1:x
        if result2008(i,1) == aws_seoul(j,1)
            result_seoul(j,6) = result2008(i,2);
        end
    end
end

[x y] = size(result2012);
for j = 1:35
    for i = 1:x
        if result2012(i,1) == aws_seoul(j,1)
            result_seoul(j,8) = result2012(i,2);
        end
    end
end

for i = 1:35
    result_seoul(i,1) = aws_seoul(i,1);
    result_seoul(i,3) = 8760;
    result_seoul(i,5) = 8760;
    result_seoul(i,7) = 8784;
    result_seoul(i,9) = 8784;
end

clear x y i j