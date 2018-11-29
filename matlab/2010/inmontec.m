function [erf_mi_rand,erf_te_rand]=inmontec(data_set,num_bin,iteration,max_tau,confidence,core)

%Information Monte Carlo Simulator
%Using randsample function, calculate mutual informaion and transfer
%entropy given random sample data set. 
%edited by Jewoo Hong(2011-04-11)
%edited by Juyeol Yun(2011-04-28)
%use parallel compute

[~,num_var]=size(data_set);
erf_mi_rand=zeros(num_var,num_var);
erf_te_rand=zeros(num_var*(max_tau),num_var);

if core < 1
    error('input argument should be positive integer')
end

matlabpool(core)

% firt random sampling

parfor i=1:iteration
    mi_sample(:,:,i)=minfo(rsam(data_set),num_bin);
    te_sample(:,:,i)=reshape(transen(rsam(data_set),num_bin,max_tau),...
                              num_var*max_tau,num_var);
end

%iteration confidence
%confidence interval 99% 
%sample error 1%
% parfor i=1:num_var
%     std_mi_sample(i)=nanstd(nanmean(mi_sample(:,i,:),3));
%     std_te_sample(i)=nanstd(nanmean(te_sample(:,i,:),3));
% end
% confidence_iteration=(2.54^2*max(max(std_mi_sample(:)),max(std_te_sample(:)))^2)/(0.01)^2;
% 
%  if iteration < confidence_iteration
%      parfor i=1:ceil(confidence_iteration-iteration)
%          add_mi_sample(:,:,i)=minfo(rsam(data_set),num_bin);
%          add_te_sample(:,:,i)=reshape(transen(rsam(data_set),num_bin,max_tau),...
%                                   num_var*max_tau,num_var);
%      end
%      mi_sample=cat(3,mi_sample,add_mi_sample);
%      te_sample=cat(3,te_sample,add_te_sample);
%  end

parfor i=1:num_var
    erf_mi_rand(:,i)=norminv(confidence,nanmean(mi_sample(:,i,:),3),nanstd(mi_sample(:,i,:),0,3));
    erf_te_rand(:,i)=norminv(confidence,nanmean(te_sample(:,i,:),3),nanstd(te_sample(:,i,:),0,3));
end

matlabpool close;

erf_te_rand=reshape(erf_te_rand,num_var,num_var,max_tau);
%csvwrite('iteration.csv',cat(1,iteration,confidence_iteration));







