 lai_co = zeros(1822,2);
 
 lai(115,2) = 0.9;
 
 for i = 1:1822
     lai_co(i,1) = i;
     if ((find(lai(:,1)==i)>0)==1)
         temp = find(lai(:,1)==i);
         lai_co(i,2) = lai(temp,2);
     else
         lai_co(i,2) = -999;
     end
 end
 clear i temp
 
 for i = 1:1822
     if lai_co(i,2) == -999
         lai_co(i,2) = lai_co(i-1,2) + (y2-y1)/(x2-x1);
         if lai_co(i,2)<0.1
             lai_co(i,2) = 0.1;
         end
     else
         if i == 1822
         else
             x1 = i;
             y1 = lai_co(i,2);
             for j = 4:8
                 if (lai_co(i+j,2) >0)
                     x2 = i+j;
                     y2 = lai_co(i+j,2);
                 end
             end
         end
     end
 end
 clear i x1 x2 y1 y2 j

 lai_co(1823,1) = 1823;
 lai_co(1824,1) = 1824;
 lai_co(1825,1) = 1825;
 lai_co(1826,1) = 1826;
 lai_co(1827,1) = 1827;
 lai_co(1823,2) = 0.2;
 lai_co(1824,2) = 0.2;
 lai_co(1825,2) = 0.2;
 lai_co(1826,2) = 0.2;
 lai_co(1827,2) = 0.2;
  
 lai_30 = zeros(1827*48,1);
 for i = 1:48
     for j = 1:1827
         lai_30((j-1)*48 + i,1) = lai_co(j,2);
     end
 end
 clear i j