%This program show process network by adjancency matrix  
%that is calculated mutual information and transfer entropy
%Data is verification by Monte-Carlo simulation.
%CSSL, SNU, 2011-04

tic;
format short

disp(...
    'type your data file full name,It should be csv file')
data_set=importdata(input('filename= ','s'));

disp(...
    'type average peroid.(unit : days)')
period=input('period_avg= ');

disp(...
    'type timestep.(unit : second)')
timestep=input('timestep= ');

disp(...
    'type average bin')
num_bin=input('num_bin= ');

disp(...
    'type maximum time lag')
max_tau=input('max_tau= ');

disp(...
    'if you want to calculate Monte-Carlo simulation, type iteration number larger than zero')
iteration=input('iteration= ');

if (iteration > 0)
    disp(...
        '< caution > This value start a paralle computing system, It should be positive integer')
    disp('If you do not want to use this value, type 1.')
    core_size=input('core_size= ');

    disp(...
        'type confidence interval you want.(0 < confidenc < 1)')
    confidence=input('confidence= ');
    if confidence>1 || confidence<0
        error('confidence should be from 0 to 1')
    end
end

if find(data_set==-99999)
    data_set=replace(data_set,-99999,NaN);
end

num_var=size(data_set,2);
te_result=zeros(num_var,num_var,max_tau+1);
erf_te_rand=zeros(num_var,num_var,max_tau+1);
final_result=zeros(max_tau+1,4*num_var+1);
w_rmi=zeros(num_var,num_var);

ano_data_set=anomaly(data_set,period,timestep);
mi_result=minfo(ano_data_set,num_bin);
te=transen(ano_data_set,num_bin,max_tau);

if (iteration > 0)
    [erf_mi_rand,erf_te]=inmontec(ano_data_set,num_bin,iteration,max_tau,confidence,core_size);
    nom_w_mi=ATz(mi_result./log10(num_bin),erf_mi_rand./log10(num_bin));
    [nom_w_te,nom_w_te_i]=ATz(te,erf_te);
     
    te=cat(3,mi_result,te);
    erf_te=cat(3,erf_mi_rand,erf_te);
    
    for l=1:max_tau+1
        te_result(:,:,l)=te(:,:,l)./mi_result;
        erf_te_rand(:,:,l)=erf_te(:,:,l)./erf_mi_rand;
    end
    
    thr_te=permute(erf_te_rand,[3,1,2]);
    ch_te_result=permute(te_result,[3,1,2]);
    
    for k=1:num_var
        w_rmi(:,k)=nom_w_mi(:,k)./nom_w_mi(k,k);
    end
 
    for i=1:num_var
        for j=1:num_var
            final_result=horzcat(permute(0:max_tau,[2,1]),ch_te_result(:,j,i),ch_te_result(:,i,j),thr_te(:,j,i),thr_te(:,i,j));
        end
        filename=strcat('final_result_',num2str(i),'.csv');
        csvwrite(filename,final_result);
    end
            
    csvwrite('ano_data_set.csv',ano_data_set);
    csvwrite('mi_result.csv',mi_result);
    csvwrite('thr_mi.csv',erf_mi_rand);
    csvwrite('nom_w_mi.csv',nom_w_mi*100);
    csvwrite('nom_w_te.csv',nom_w_te./mi_result);
    csvwrite('nom_w_te_index.csv',nom_w_te_i);
    csvwrite('w_rmi.csv',w_rmi*100);
end

clear i j l k core_size filename final_result thr_te ch_te_result
toc;


