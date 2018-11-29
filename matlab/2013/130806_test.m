%% 130806, for test
[n_st n_year n_data n_var] = size(raw_aws_seoul);

raw_seoul = zeros(n_st,n_year*n_data,n_var);
for i = 1:n_st
    for j = 1:n_year
        for k = 1:n_data
            for l = 1:n_var
                raw_seoul(i,(j-1)*n_data+k,l) = raw_aws_seoul(i,j,k,l);
            end
        end
    end
end
clear i j k l


%% 0. Gap (replaced by -999)
result_gap = zeros(n_st,(n_var-1));
for i = 1:n_st
    for j = 1:(n_var-1)
        temp = find(raw_seoul(i,:,j+1)==-999);
        [x, y] = size(temp);
        if x==0
            result_gap(i,j) =0;
        elseif y==0
            result_gap(i,j) =0;
        else
            if x>y
                result_gap(i,j) =x;
            elseif x==y
                result_gap(i,j) =x;
            else
                result_gap(i,j) =y;
            end
        end
    end
end
for i = 1:n_st
    for j = 1:(n_var-1)
        result_gap(i,j) = result_gap(i,j)-48; %remove 29th Feb. in 2006, 2007
        if i ==25
            result_gap(i,j) = result_gap(i,j)- 8760-8760-8784; % remove 2006, 2007, 2008(no data)
        elseif i ==28
            result_gap(i,j) = result_gap(i,j)- 8760-8760-8784; % remove 2006, 2007, 2008(no data)
        elseif i ==31
            result_gap(i,j) = result_gap(i,j)- 8760-8760-8784; % remove 2006, 2007, 2008(no data)
        end
    end
end
clear i j x y temp 


%% 1. Physical Range Test (replaced by -998)
mod_seoul = raw_seoul;
result_range = zeros(n_st,(n_var-1)*2); % min and max
range_var = [-40, 100; 0, 150; 0, 100; 800, 1040; 0, 70; 0 360]; 
% standard for range test
% references are according to "Weather Observation Standardization Act"

for i = 1:n_st
    for j = 1:(n_var-1)
        for k = 1:2
            if k ==1
                temp = find(raw_seoul(i,:,j+1)<range_var(j,k));
            else
                temp = find(raw_seoul(i,:,j+1)>range_var(j,k));
            end
            [x, y] = size(temp);
            if x==0
                result_range(i,(j-1)*2+k) =0;
            elseif y==0
                result_range(i,(j-1)*2+k) =0;
            else
                if x>y
                    for l = 1:x
                        if mod_seoul(i,temp(l,1),j+1) ~= -999
                            mod_seoul(i,temp(l,1),j+1) = -998;
                        end
                    end
                    result_range(i,(j-1)*2+k) =x;
                elseif x==y
                    result_range(i,(j-1)*2+k) =x;
                    if mod_seoul(i,temp(x,1),j+1) ~= -999
                        mod_seoul(i,temp(x,1),j+1) = -998;
                    end
                else
                    for l = 1:y
                        if mod_seoul(i,temp(1,l),j+1) == -999
                        else
                            mod_seoul(i,temp(1,l),j+1) = -998;
                        end
                    end
                    result_range(i,(j-1)*2+k) =y;
                end
            end
        end
    end
end

% remove the number of the gap
for i = 1:n_st
    for j = 1:(n_var-1)
        for k = 1:1
            result_range(i,(j-1)*2+k) = result_range(i,(j-1)*2+k)-result_gap(i,j)-48;
        end
        if i ==25
            result_range(i,(j-1)*2+k) = result_range(i,(j-1)*2+k)-8760-8760-8784;
        elseif i ==28
            result_range(i,(j-1)*2+k) = result_range(i,(j-1)*2+k)-8760-8760-8784;
        elseif i ==31
            result_range(i,(j-1)*2+k) = result_range(i,(j-1)*2+k)-8760-8760-8784;
        end
    end
end
clear i j k l temp x y range_var


%% 2. Internal Consistency Check (replaced by -997)

result_consistency = zeros(n_st,2);

for i = 1:n_st
    j=5; % wind-speed
    temp = find(mod_seoul(i,:,j+1)==0);
    [x, y] = size(temp);
    if y>=1
        for k=1:y
            if ((raw_seoul(i,temp(1,k),j+2)~=0)&&(raw_seoul(i,temp(1,k),j+2)~=-999))&&(raw_seoul(i,temp(1,k),j+2)~=-998)
                result_consistency(i,1) = result_consistency(i,1) +1;
                mod_seoul(i,temp(1,k),j+1) = -997;
                mod_seoul(i,temp(1,k),j+2) = -997;
            end
        end
    end
    
    temp = find(mod_seoul(i,:,j+2)==0);
    [x, y] = size(temp);
    if y>=1
        for k=1:y
            if ((raw_seoul(i,temp(1,k),j+1)~=0)&&(raw_seoul(i,temp(1,k),j+1)~=-999))&&(raw_seoul(i,temp(1,k),j+1)~=-998)
                result_consistency(i,2) = result_consistency(i,2) +1;
                mod_seoul(i,temp(1,k),j+1) = -997;
                mod_seoul(i,temp(1,k),j+2) = -997;
            end
        end
    end
    
end
clear i j k x y temp


%% 3. Step Test (replaced by -996)

result_step = zeros(n_st, 4);
temp_std = zeros(n_st,4);
temp_n = zeros(n_st,4);
for i = 1:n_st
    temp = 0;
    temp1 = 0;
    for j = 1:(n_data*4-1)
        k=1;
        if mod_seoul(i,j,k+1)>-990
            temp = temp + 1;
            temp1(temp,1) = mod_seoul(i,j,k+1);
        end
    end
    temp_std(i,1) = std(temp1);
    temp_n(i,1) = temp;
    
    temp = 0;
    temp1 = 0;
    for j = 1:(n_data*4-1)
        k=3;
        if mod_seoul(i,j,k+1)>-990
            temp = temp + 1;
            temp1(temp,1) = mod_seoul(i,j,k+1);
        end
    end
    temp_std(i,2) = std(temp1);
    temp_n(i,2) = temp;    
    
    temp = 0;
    temp1 = 0;
    for j = 1:(n_data*4-1)
        k=4;
        if mod_seoul(i,j,k+1)>-990
            temp = temp + 1;
            temp1(temp,1) = mod_seoul(i,j,k+1);
        end
    end
    temp_std(i,3) = std(temp1);
    temp_n(i,3) = temp;    
    
    temp = 0;
    temp1 = 0;
    for j = 1:(n_data*4-1)
        k=5;
        if mod_seoul(i,j,k+1)>-990
            temp = temp + 1;
            temp1(temp,1) = mod_seoul(i,j,k+1);
        end
    end
    temp_std(i,4) = std(temp1);        
    temp_n(i,4) = temp;
    
end

for i = 1:n_st
    temp = 0;
    temp1 = 0;
    for j = 2:(n_data*4-1)
        k=1;
        if mod_seoul(i,j,k+1)>-990
            if (mod_seoul(i,j-1,k+1)>-990)&&(mod_seoul(i,j+1,k+1)>-990)
                temp1 = abs(mod_seoul(i,j-1,k+1)-mod_seoul(i,j,k+1))+abs(mod_seoul(i,j+1,k+1)-mod_seoul(i,j,k+1));
                if temp1>4*temp_std(i,1)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end
            elseif mod_seoul(i,j-1,k+1)>-990
                temp1 = abs(mod_seoul(i,j-1,k+1)-mod_seoul(i,j,k+1));
                if temp1>2*temp_std(i,1)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end
            elseif mod_seoul(i,j+1,k+1)>-990
                temp1 = abs(mod_seoul(i,j+1,k+1)-mod_seoul(i,j,k+1));
                if temp1>2*temp_std(i,1)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end                
            end
        end
    end
    result_step(i,1) = temp;
    
    temp = 0;
    temp1 = 0;
    for j = 2:(n_data*4-1)
        k=3;
        if mod_seoul(i,j,k+1)>-990
            if (mod_seoul(i,j-1,k+1)>-990)&&(mod_seoul(i,j+1,k+1)>-990)
                temp1 = abs(mod_seoul(i,j-1,k+1)-mod_seoul(i,j,k+1))+abs(mod_seoul(i,j+1,k+1)-mod_seoul(i,j,k+1));
                if temp1>4*temp_std(i,2)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end
            elseif mod_seoul(i,j-1,k+1)>-990
                temp1 = abs(mod_seoul(i,j-1,k+1)-mod_seoul(i,j,k+1));
                if temp1>2*temp_std(i,2)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end
            elseif mod_seoul(i,j+1,k+1)>-990
                temp1 = abs(mod_seoul(i,j+1,k+1)-mod_seoul(i,j,k+1));
                if temp1>2*temp_std(i,2)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end                
            end
        end        
    end
    result_step(i,2) = temp;    
    
    temp = 0;
    temp1 = 0;
    for j = 2:(n_data*4-1)
        k=4;
        if mod_seoul(i,j,k+1)>-990
            if (mod_seoul(i,j-1,k+1)>-990)&&(mod_seoul(i,j+1,k+1)>-990)
                temp1 = abs(mod_seoul(i,j-1,k+1)-mod_seoul(i,j,k+1))+abs(mod_seoul(i,j+1,k+1)-mod_seoul(i,j,k+1));
                if temp1>4*temp_std(i,3)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end
            elseif mod_seoul(i,j-1,k+1)>-990
                temp1 = abs(mod_seoul(i,j-1,k+1)-mod_seoul(i,j,k+1));
                if temp1>2*temp_std(i,3)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end
            elseif mod_seoul(i,j+1,k+1)>-990
                temp1 = abs(mod_seoul(i,j+1,k+1)-mod_seoul(i,j,k+1));
                if temp1>2*temp_std(i,3)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end                
            end
        end
    end
    result_step(i,3) = temp;    
    
    temp = 0;
    temp1 = 0;
    for j = 2:(n_data*4-1)
        k=5;
        if mod_seoul(i,j,k+1)>-990
            if (mod_seoul(i,j-1,k+1)>-990)&&(mod_seoul(i,j+1,k+1)>-990)
                temp1 = abs(mod_seoul(i,j-1,k+1)-mod_seoul(i,j,k+1))+abs(mod_seoul(i,j+1,k+1)-mod_seoul(i,j,k+1));
                if temp1>4*temp_std(i,4)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end
            elseif mod_seoul(i,j-1,k+1)>-990
                temp1 = abs(mod_seoul(i,j-1,k+1)-mod_seoul(i,j,k+1));
                if temp1>2*temp_std(i,4)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end
            elseif mod_seoul(i,j+1,k+1)>-990
                temp1 = abs(mod_seoul(i,j+1,k+1)-mod_seoul(i,j,k+1));
                if temp1>2*temp_std(i,4)
                    temp = temp +1;
                    mod_seoul(i,j,k+1) = -996;
                end                
            end
        end
    end
    result_step(i,4) = temp;
end
clear i j k temp1 temp temp_n temp_std 



