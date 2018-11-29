function z=anomaly(data_set, num_var, num_data, period_avg, timestep)
% anomaly dataset
% data_set : input time-series data set
% period_avg : default : 5days
% timestep : default : 1800secs

% data_set information
num_data_day=24*3600/timestep;
num_day=num_data(1,1)/num_data_day;
avg_sample=0;
ano_data_set=zeros(num_data-num_data_day*(period_avg-1),num_var);

% anomaly dataset
for i=1:num_var
    for j=(period_avg-1)*num_data_day+1:num_data
        if data_set(j,i)~=-99999
            z=0;
            for k=1:period_avg
                if data_set(j-(k-1)*num_data_day,i)==-99999
                else
                    avg_sample=avg_sample+data_set(j-(k-1)*num_data_day,i);
                    z=z+1;
                end
            end
            avg_sample=avg_sample/z;
            ano_data_set(j-(period_avg-1)*num_data_day,i)=data_set(j,i)-avg_sample;
            avg_sample=0;
        else 
            ano_data_set(j-(period_avg-1)*num_data_day,i)=-99999;
        end
    end
end

% return
z=ano_data_set;

