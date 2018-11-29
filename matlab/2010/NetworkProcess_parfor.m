tic
%option 
num_data_day=24;
period_avg=7;
timestep=3600;
num_bin=32;
max_tau=24;
cycle=1000;
% data_set information
size_data=size(data_set);
num_var=size_data(1,2);
num_data=size_data(1,1);
% clear size_data

tic
% anomaly
ano_data_set=anomaly(data_set, num_var, num_data, period_avg, timestep);
clear period_avg timestep num_data
% anomaly information
size_data=size(ano_data_set);
num_ano_data=size_data(1,1);
num_day=num_ano_data(1,1)/num_data_day;
% clear size_data num_data_day
toc

% tic
matlabpool(4);
% toc

% 130426_ anomaly cancel for electricity experiment
% ano_data_set = data_set;

tic
% shannon entropy
H=zeros(num_var,1);
for i=1:num_var
    H(i)=shannonentropy(ano_data_set, i, num_bin, num_ano_data);
end
toc

tic
% mutual information
mi_result=zeros(num_var,num_var);
for i=1:num_var
    for j=1:num_var
        mi_result(i,j)=mutualinformation(ano_data_set, i, j, num_bin, num_ano_data);
    end
end
toc

tic
% transfer entropy
h=waitbar(0,'Transfer Entropy, Please Wait...');
te_result=zeros(num_var,num_var,max_tau);
for i=1:num_var
    waitbar(i/num_var)
    for j=1:num_var
        parfor tau=1:max_tau
            te_result(i,j,tau)=transferentropy(ano_data_set, i, j, tau, num_bin, num_ano_data);
        end
    end
end
close(h);
toc

tic
% monte carlo simulation
h=waitbar(0,'MI. Monte Carlo, Please Wait...');
mi_mean_rand=zeros(num_var,num_var); %mean
mi_sigma_rand=zeros(num_var,num_var); %std
mi_sample_rand=zeros(num_var,num_var,cycle); %result
for i=1:num_var
    waitbar(i/num_var)
    for j=1:num_var
        parfor c=1:cycle
            mi_sample_rand(i,j,c)=montemutual(ano_data_set, i, j, num_bin, num_ano_data);
        end
    end
end
close(h);
toc
for i=1:num_var
    for j=1:num_var
        mi_mean_rand(i,j)=mean(mi_sample_rand(i,j,:));
        mi_sigma_rand(i,j)=std(mi_sample_rand(i,j,:));
    end
end
clear mi_sample_rand
tic
h=waitbar(0,'TE. Monte Carlo, Please Wait...');
te_mean_rand=zeros(num_var,num_var,max_tau); %mean
te_sigma_rand=zeros(num_var,num_var,max_tau); %std
te_sample_rand=zeros(num_var,num_var,max_tau,cycle); %result
for i=1:num_var
    waitbar(i/num_var)
    for j=1:num_var
        for tau=1:max_tau
            parfor c=1:cycle
                te_sample_rand(i,j,tau,c)=montetrans(ano_data_set, i, j, tau, num_bin, num_ano_data);
            end
        end
    end
end
close(h);
toc
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau
            te_mean_rand(i,j,k)=mean(te_sample_rand(i,j,k,:));
            te_sigma_rand(i,j,k)=std(te_sample_rand(i,j,k,:));
        end
    end
end
clear te_sample_rand


matlabpool close;
toc