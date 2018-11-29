l = 98;

temp = 0;
temp_l = 0;
for i = 1:l
    temp = size(total_result{i,2});
    temp_l = temp_l + temp(1,1);
end 
clear i
final_result = zeros(temp_l,54);
temp = 0;
temp_l = 0;
for i = 1:l
    temp = size(total_result{i,2});
    for j = (temp_l+1):(temp_l+temp(1,1))
        final_result(j,:) = total_result{i,2}(j-temp_l,:);
    end
    clear j
    temp_l = temp_l + temp(1,1);
end
clear i temp temp_l l 