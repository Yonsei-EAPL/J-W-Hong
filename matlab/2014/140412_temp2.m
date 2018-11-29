%%
a=sort(output,'descend');

for i=1:4001
   b(i,:)=sort(a(i,:),'descend');
end
d=sort(c,'descend');

%%
[size_x size_y]=size(output);
temp_data=0;
a=1;
b=1;
c=1;
d=1;

for i=1:size_x
    for j=1:size_y
        if output(i,j)>10^(-1)
        temp_data(a,1)=output(i,j);
        a=a+1;
        elseif output(i,j)>10^(-2)
        temp_data(b,2)=output(i,j);
        b=b+1;
        elseif output(i,j)>10^(-3)
        temp_data(c,3)=output(i,j);
        c=c+1;
        elseif output(i,j)>10^(-4)
        temp_data(d,4)=output(i,j);
        d=d+1;     
        end              
    end    
end

clear i j a b c d  size_x size_y

%%
temp_all = 0;
for i = 1:4001
    for j = 1 :4001
        temp_all=temp_all+output(i,j);
    end
end

temp_4number = 0;
for i = 1:4001
    for j = 1 :4001
        if (output(i,j)>10^(-4)&&(output(i,j)<10^(-3)))
            temp_4number=temp_4number+1;
        end
    end
end
clear i j

temp_5number = 0;
for i = 1:4001
    for j = 1 :4001
        if (output(i,j)>10^(-5)&&(output(i,j)<10^(-4)))
            temp_5number=temp_5number+1;
        end
    end
end
clear i j

temp_6number = 0;
for i = 1:4001
    for j = 1 :4001
        if (output(i,j)>10^(-6)&&(output(i,j)<10^(-5)))
            temp_6number=temp_6number+1;
        end
    end
end
clear i j

temp_7number = 0;
for i = 1:4001
    for j = 1 :4001
        if (output(i,j)>10^(-7)&&(output(i,j)<10^(-6)))
            temp_7number=temp_7number+1;
        end
    end
end
clear i j

temp_4 = 0;
for i = 1:4001
    for j = 1 :4001
        if output(i,j)>10^(-4)
            temp_4=temp_4+output(i,j);
        end
    end
end
clear i j

temp_5 = 0;
for i = 1:4001
    for j = 1 :4001
        if output(i,j)>10^(-5)
            temp_5=temp_5+output(i,j);
        end
    end
end
clear i j

temp_6 = 0;
for i = 1:4001
    for j = 1 :4001
        if output(i,j)>10^(-6)
            temp_6=temp_6+output(i,j);
        end
    end
end
clear i j

temp_7 = 0;
for i = 1:4001
    for j = 1 :4001
        if output(i,j)>10^(-7)
            temp_7=temp_7+output(i,j);
        end
    end
end
clear i j


