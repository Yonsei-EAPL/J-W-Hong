    matlabpool(4);
        
    %option 
    num_data_day=24;
    period_avg=7;
    timestep=3600;
    max_tau=24;
    cycle=1000;
    % data_set information
    size_data=size(data_set);
    num_var=size_data(1,2);
    num_data=size_data(1,1);
    % clear size_data
    
    max_bin = 1500;
    H=zeros(max_bin, num_var);
    mi_result=zeros(max_bin, num_var, num_var);
for nbin=1:max_bin

    num_bin=nbin;

    % anomaly
    ano_data_set=anomaly(data_set, num_var, num_data, period_avg, timestep);
%     clear period_avg timestep num_data
    % anomaly information
    size_data=size(ano_data_set);
    num_ano_data=size_data(1,1);
    num_day=num_ano_data(1,1)/num_data_day;
    % clear size_data num_data_day

    % 130426_ anomaly cancel for electricity experiment
    % ano_data_set = data_set;

    % shannon entropy
    for i=1:num_var
        H(nbin, i)=shannonentropy(ano_data_set, i, num_bin, num_ano_data);
    end

    % mutual information
    for i=1:num_var
        for j=1:num_var
            mi_result(nbin, i,j)=mutualinformation(ano_data_set, i, j, num_bin, num_ano_data);
        end
    end
 
end

matlabpool close;

temp = zeros(max_bin,1);
temp2 = zeros(max_bin,1);
for i = 1:max_bin
    temp(i,1) = mi_result(i,1,1);
    temp2(i,1) = mi_result(i,2,1);
end
plot(temp)
hold on 
plot(temp2,'r')
hold on
for i = 1:max_bin
    temp2(i,1) = temp2(i,1)/temp(i,1);
end
plot(temp2,'c')
