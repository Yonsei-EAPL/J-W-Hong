n=0;
temp2=0;
tempd=0;

for i = 1:5875
    if i ==1
        p=1;
        n=1;
        temp2=a(1,3);
        tempd=a(1,4);
    else
        if (temp2==a(i,3))&&(n<5)
            n=n+1;
            if tempd>a(i,4)
                tempd = a(i,4);
            end
        elseif n>3
            result(p,1) = a(i,1);
            result(p,2) = a(i,2);
            result(p,3) = temp2;
            result(p,4) = tempd;
            result(p,5) = a(i,5);
            
            p=p+1;
            n=1;
            temp2 = a(i,3);
            tempd = a(i,4);
        end
    end
end




