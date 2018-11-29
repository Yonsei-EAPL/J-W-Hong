B = zeros(2501*2501,1);
for i = 1:25020004/4
	for j = 1:4
        if j ==1
            temp_sign =0;
            temp_exp = 0;
            temp_res = 1;
            if A((i-1)*4+j,1)>128
                temp_sign = 1;
                temp_exp = (A((i-1)*4+j,1)-128)*2;
            else
                temp_exp = (A((i-1)*4+j,1))*2;
            end
        end
        if j ==2
            if A((i-1)*4+j,1)>128
                temp_exp = temp_exp + 1;
                temp_exp = temp_exp - 127;
                temp_res = temp_res + (A((i-1)*4+j,1)-128)*2^(-8);
            else
                temp_exp = temp_exp - 127;
                temp_res = temp_res + (A((i-1)*4+j,1)-128)*2^(-8);
            end
        end
        if j ==3
            temp_res = temp_res + A((i-1)*4+j,1)*2^(-16);
        end
        if j ==4
            temp_res = temp_res + A((i-1)*4+j,1)*2^(-24);
            if temp_sign ==1
                B(i,1) = (-1)*2^(temp_exp)*temp_res;
            else
                B(i,1) = 2^(temp_exp)*temp_res;
            end
        end
    end
    
end
clear i j temp_sign temp_exp temp_res

C = zeros(2501,2501);
for i = 1:2501
    for j = 1:2501
        C(i,j) = B((i-1)*2501+j,1);
    end
end
clear i j 
C_inv = C^(-1);

O = eye(2501);
for i = 1:2501
    for j = 1:2501
        if O(i,j) >0
            O(i,j) = 51^2;
        end
    end
end
clear i j
Result2 = O^(-1);

Result1 = (C_inv + Result2)^(-1);

temp = 0;
for i = 1:2501
    for j = 1:2501
        if i ==j
            temp = temp + C(i,j);
        end
    end
end
clear i j temp