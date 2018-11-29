final1 = zeros(35,3); % mean, standard error, number
final2 = zeros(35,3);
final3 = zeros(35,3);
final4 = zeros(35,3);
final5 = zeros(35,3);

for i = 1:35
    for j = 87649:113952
%         if (time(j,4)==1)&&(time(j,5)==1)
        if time(j,5)==1
            if result1(j,i)>0
                final1(i,3) = final1(i,3)+1;
            end
            if result2(j,i)>0
                final2(i,3) = final2(i,3)+1;
            end
            if result3(j,i)>0
                final3(i,3) = final3(i,3)+1;
            end
            if result4(j,i)>0
                final4(i,3) = final4(i,3)+1;
            end
            if result5(j,i)>0
                final5(i,3) = final5(i,3)+1;
            end            
        end
    end
end
clear i j 

for i = 1:35
    temp1 = zeros(final1(i,3),1);
    temp2 = zeros(final2(i,3),1);
    temp3 = zeros(final3(i,3),1);
    temp4 = zeros(final4(i,3),1);
    temp5 = zeros(final5(i,3),1);
    temp11 = 0;
    temp22 = 0;
    temp33 = 0;
    temp44 = 0;
    temp55 = 0;

    for j = 87649:113952
%         if (time(j,4)==1)&&(time(j,5)==1)
        if time(j,5)==1
            if result1(j,i)>0
                temp11 = temp11+1;
                temp1(temp11,1) = result1(j,i);
            end
            if result2(j,i)>0
                temp22 = temp22+1;
                temp2(temp22,1) = result2(j,i);
            end
            if result3(j,i)>0
                temp33 = temp33+1;
                temp3(temp33,1) = result3(j,i);
            end
            if result4(j,i)>0
                temp44 = temp44+1;
                temp4(temp44,1) = result4(j,i);
            end
            if result5(j,i)>0
                temp55 = temp55+1;
                temp5(temp55,1) = result5(j,i);
            end            
        end
    end
    
    final1(i,1) = mean(temp1);
    final1(i,2) = std(temp1)/(final1(i,3))^(0.5);
    final2(i,1) = mean(temp2);
    final2(i,2) = std(temp2)/(final2(i,3))^(0.5);    
    final3(i,1) = mean(temp3);
    final3(i,2) = std(temp3)/(final3(i,3))^(0.5);    
    final4(i,1) = mean(temp4);
    final4(i,2) = std(temp4)/(final4(i,3))^(0.5);    
    final5(i,1) = mean(temp5);
    final5(i,2) = std(temp5)/(final5(i,3))^(0.5);        
end
clear i j 