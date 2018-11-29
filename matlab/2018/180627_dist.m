output= zeros(230,5);
output(:,1)=input(:,1);
for i = 1:230
    
    dist = zeros(569,2);
    for j = 1:569
        dist(j,1) = ((data(j,2)-input(i,2))^2 + (data(j,3)-input(i,3))^2)^0.5;
        dist(j,2) = data(j,1);
    end
    clear j
    dist = sortrows(dist,1);
    
    for j = 1:4
        output(i,j+1)= dist(j,2);
    end
    clear j
end
clear i 