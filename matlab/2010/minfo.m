function mi_result=minfo(data_set,num_bin)

%This fuction calculate mutual informaiton given data_set.
%The result is shown 
%Number_bin is binary representation of stored integer of object
%declare variables
%edited by Jewoo Hong(20num_bin-04-11)
%edited by Juyeol Yun(2011-04-28)
[~,num_var]=size(data_set);
pr_data=zeros(num_bin,num_var);
mi_result=zeros(num_var,num_var);


for i=1:num_var
    pr_data(:,i)=hist(data_set(:,i),num_bin)/sum(hist(data_set(:,i),num_bin));
end

for j=1:num_var
    for k=1:num_var
        jpr_data=hist3(horzcat(data_set(:,j),data_set(:,k)),[num_bin,num_bin])...
                 /nansum(nansum(hist3(horzcat(data_set(:,j),data_set(:,k)),[num_bin,num_bin])));       
        mi_result(j,k)=(-1)*nansum(pr_data(:,j).*log10(pr_data(:,j)))+(-1)*nansum(pr_data(:,k).*log10(pr_data(:,k)))-...
                       (-1)*nansum(nansum(jpr_data.*log10(jpr_data)));
    end
end





