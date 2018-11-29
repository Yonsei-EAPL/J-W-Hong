function [Weighted_Adjacency_data_set,Weighted_Adjacenccy_position]=ATz(data_set,erf_data_set)

%This function make a network adjacency matrix.
%A network adjacency matrix is weighted cut matrix.
% 2011. 3. 28 edited by Juyeol 
num_var=size(data_set,2);
data_set( data_set < erf_data_set & abs(data_set-erf_data_set) > 1e-5)=0;

[Weighted_Adjacency_data_set,position] = nanmax(data_set,[],3);
Weighted_Adjacenccy_position=nan(num_var*num_var,4);

k=1;
for i=1:15
    for j=1:15
        if isempty(find(data_set(i,j,:)>0,1,'first')/2)==0
            Weighted_Adjacenccy_position(k,:)=[find(data_set(i,j,:)>0,1,'first')/2,...
                find(data_set(i,j,:)>0,1,'last')/2,length(find(data_set(i,j,:)>0)),position(i,j)/2];
        end
        k=k+1;
    end
end







    





        
