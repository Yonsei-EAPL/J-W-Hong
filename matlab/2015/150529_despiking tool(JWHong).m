des_sigma = 3.5; % thresohold standard deviation value for despiking
des_window = 60; % window width in minute
des_window = des_window * 60 * 10;  

% despiking
des_n = 0;
for j = 1:(length(co2)-des_window+1)
    n=0;
    des_n_temp = 1;
    des_n_before = 0;
    n_start = j;
    n_end = j+des_window-1;    
    while des_n_temp>des_n_before
        des_n_before = des_n_temp;
        n=n+1;
        des_std = abs(std(co2(n_start:n_end,2)));
        des_mean = mean(co2(n_start:n_end,2));
        for k = 1:des_window
            if abs(co2(j+k-1,2)-des_mean)>(des_sigma+(n-1)*0.1)*des_std
                if k==1
                    co2(j+k-1,2)=mean(co2(n_start:n_end,2));
                    co2(j+k-1,3)=co2(j+k-1,3)+1;
                    des_n = des_n +1;
                    des_n_temp = des_n_temp +1;                    
                elseif k==des_window
                    co2(j+k-1,2)=mean(co2(n_start:n_end,2));
                    co2(j+k-1,3)=co2(j+k,3)+1;
                    des_n = des_n +1;
                    des_n_temp = des_n_temp +1;                    
                else
                    co2(j+k-1,2)=mean(co2(n_start:n_end,2));
                    co2(j+k-1,3)=co2(j+k,3)+1;
                    des_n = des_n +1;
                    des_n_temp = des_n_temp +1;
                end
            end
        end
    end
end
%clear des_sigma des_window des_n des_n_temp des_n_before
%clear j k n
