function rand_data_set=rsam(data_set)

%This function make a random data set using randsample function
%parallel compute


[num_data,num_var]=size(data_set);
rand_data_set=zeros(num_data,num_var);

for i=1:num_var
    rdata_set=data_set(:,i);
    rand_data_set(:,i)=randsample(rdata_set,num_data);
end
