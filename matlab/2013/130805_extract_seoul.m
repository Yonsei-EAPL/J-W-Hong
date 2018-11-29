%%
raw_aws_seoul = zeros(34,4,8784,7);
for i = 1:34
    for j = 1:4
        for k = 1:8784
            for l = 2:7
                raw_aws_seoul(i,j,k,l) = -999;
            end
            raw_aws_seoul(i,j,k,1) = time(k,1);
        end
    end
end
clear i j k l

%% for 2006
clear result
n_year = 8760;
[x y] = size(aws2006);
for i = 1:x
    aws2006(i,1) = mod(aws2006(i,1),1000000);
end
clear i

temp = zeros(34,1);
for i = 1:x
    for j = 1:34
        if aws_seoul(j,1) == aws2006(i,2)
            temp(j,1) = temp(j,1) + 1;
            temp1 = find(time==aws2006(i,1));
            raw_aws_seoul(j,1,temp1,2) = aws2006(i,3);
            raw_aws_seoul(j,1,temp1,3) = aws2006(i,5);
            raw_aws_seoul(j,1,temp1,4) = aws2006(i,6);
            raw_aws_seoul(j,1,temp1,5) = aws2006(i,7);
            raw_aws_seoul(j,1,temp1,6) = aws2006(i,9);
            raw_aws_seoul(j,1,temp1,7) = aws2006(i,10);            
        end
    end
end
clear i j 
clear temp1
clear x y 
% clear aws2006

%% for 2007
clear result
n_year = 8760;
[x y] = size(aws2007);
for i = 1:x
    aws2007(i,1) = mod(aws2007(i,1),1000000);
end
clear i

temp = zeros(34,1);
for i = 1:x
    for j = 1:34
        if aws_seoul(j,1) == aws2007(i,2)
            temp(j,1) = temp(j,1) + 1;
            temp1 = find(time==aws2007(i,1));
            raw_aws_seoul(j,2,temp1,2) = aws2007(i,3);
            raw_aws_seoul(j,2,temp1,3) = aws2007(i,5);
            raw_aws_seoul(j,2,temp1,4) = aws2007(i,6);
            raw_aws_seoul(j,2,temp1,5) = aws2007(i,7);
            raw_aws_seoul(j,2,temp1,6) = aws2007(i,9);
            raw_aws_seoul(j,2,temp1,7) = aws2007(i,10);            
        end
    end
end
clear i j 
clear temp1
clear x y aws2007

%% for 2008
clear result
n_year = 8784;
[x y] = size(aws2008);
for i = 1:x
    aws2008(i,1) = mod(aws2008(i,1),1000000);
end
clear i

temp = zeros(34,1);
for i = 1:x
    for j = 1:34
        if aws_seoul(j,1) == aws2008(i,2)
            temp(j,1) = temp(j,1) + 1;
            temp1 = find(time==aws2008(i,1));
            raw_aws_seoul(j,3,temp1,2) = aws2008(i,3);
            raw_aws_seoul(j,3,temp1,3) = aws2008(i,5);
            raw_aws_seoul(j,3,temp1,4) = aws2008(i,6);
            raw_aws_seoul(j,3,temp1,5) = aws2008(i,7);
            raw_aws_seoul(j,3,temp1,6) = aws2008(i,9);
            raw_aws_seoul(j,3,temp1,7) = aws2008(i,10);            
        end
    end
end
clear i j 
clear temp1
clear x y aws2008

%% for 2012
clear result
n_year = 8784;
[x y] = size(aws2012);
for i = 1:x
    aws2012(i,1) = mod(aws2012(i,1),1000000);
end
clear i

temp = zeros(34,1);
for i = 1:x
    for j = 1:34
        if aws_seoul(j,1) == aws2012(i,2)
            temp(j,1) = temp(j,1) + 1;
            temp1 = find(time==aws2012(i,1));
            raw_aws_seoul(j,4,temp1,2) = aws2012(i,3);
            raw_aws_seoul(j,4,temp1,3) = aws2012(i,5);
            raw_aws_seoul(j,4,temp1,4) = aws2012(i,6);
            raw_aws_seoul(j,4,temp1,5) = aws2012(i,7);
            raw_aws_seoul(j,4,temp1,6) = aws2012(i,9);
            raw_aws_seoul(j,4,temp1,7) = aws2012(i,10);            
        end
    end
end
clear i j 
clear temp1
clear x y aws2012


