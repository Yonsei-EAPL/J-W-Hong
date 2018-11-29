% for SF nocturnal NEE analysis
% data 
%   1: Tk
%   2: SWC (%)
%   3: Fc

min_T = min(data(:,1));
max_T = max(data(:,1));
min_S = min(data(:,2));
max_S = max(data(:,2));
r_S = max_S - min_S;
r_T = max_T - min_T;

n = 1255;
result = zeros(10,10);
temp_sum = zeros(10,10); % SWC, Tk
temp_n = zeros(10,10); % SWC, Tk

temp_S =0;
temp_T =0;
for i = 1:n
    % find SWC address
    temp_S = ((data(i,2) - min_S) - mod((data(i,2) - min_S),(r_S/10)))/(r_S/10);
    
    temp_S = temp_S - mod(temp_S,1);
    temp_S = temp_S + 1;
    if temp_S >10
        temp_S =10;
    elseif temp_S ==0
        temp_S = 1;
    end
    % find Tk address
    temp_T = ((data(i,1) - min_T) - mod((data(i,1) - min_T),(r_T/10)))/(r_T/10);
    temp_T = temp_T - mod(temp_T,1);
    temp_T = temp_T + 1;
    if temp_T >10
        temp_T =10;
    elseif temp_T ==0
        temp_T = 1;
    end
    
    % sum
    temp_sum(temp_S, temp_T) = temp_sum(temp_S, temp_T) + data(i,3);
    temp_n(temp_S, temp_T) = temp_n(temp_S, temp_T) + 1;
end
clear i
for i = 1:10
    for j = 1:10
        result(i,j) = temp_sum(i,j)/temp_n(i,j);
    end
end
clear i j 


