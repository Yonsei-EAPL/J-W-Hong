temp=0;
temp1=0;
for i = 1:3977018
    if aws2012(i,2)==400
        temp= temp+1;
    elseif aws2012(i,2)==416
        temp1=temp1+1;
    end
end
temp
temp1
clear i

aws2012_400 = zeros(temp,12);
aws2012_416 = zeros(temp1,12);

temp=0;
temp1=0;
for i = 1:3977018
    if aws2012(i,2)==400
        temp= temp+1;
        aws2012_400(temp,:)=aws2012(i,:);
    elseif aws2012(i,2)==416
        temp1=temp1+1;
        aws2012_416(temp1,:)=aws2012(i,:);
    end    
end
clear i
