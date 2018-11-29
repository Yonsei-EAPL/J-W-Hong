result = zeros(11,6);
length = 136800;
temp = 0;

for i = 2:length
    for j = 1:11
        if x(i-1,j)==x(i,j)
            result(j,1) = result(j,1)+1;
        end
        if i>2
            if (x(i-2,j)==x(i,j))&&(x(i-1,j)==x(i,j))
                result(j,2) = result(j,2)+1;
            end
        end
        if i>3
            if ((x(i-3,j)==x(i,j))&&(x(i-2,j)==x(i,j)))&&(x(i-1,j)==x(i,j))
                result(j,3) = result(j,3)+1;
            end            
        end
        if i>4
            if (((x(i-4,j)==x(i,j))&&(x(i-3,j)==x(i,j)))&&(x(i-2,j)==x(i,j)))&&(x(i-1,j)==x(i,j))
                result(j,4) = result(j,4)+1;
            end                        
        end
        if i>5
            if ((((x(i-5,j)==x(i,j))&&(x(i-4,j)==x(i,j)))&&(x(i-3,j)==x(i,j)))&&(x(i-2,j)==x(i,j)))&&(x(i-1,j)==x(i,j))
                result(j,5) = result(j,5)+1;
            end                                    
        end
        if i>6
            if (((((x(i-6,j)==x(i,j))&&(x(i-5,j)==x(i,j)))&&(x(i-4,j)==x(i,j)))&&(x(i-3,j)==x(i,j)))&&(x(i-2,j)==x(i,j))&&(x(i-1,j)==x(i,j)))
                result(j,6) = result(j,6)+1;
                temp = temp + 1;
                add(temp,1) = i;
                add(temp,2) = j;
            end                                                
        end
    end
end
clear i j