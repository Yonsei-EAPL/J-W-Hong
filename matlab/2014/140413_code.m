%%
sort1_output=sort(output,'descend');
for i=1:4001
   sort2_output(i,:)=sort(sort1_output(i,:),'descend');
end
clear sort1_output i

%%
for i=1:2000
    for j=1:2000
        temp_output(i,j)=sort2_output(i,j);
    end
end
clear i j
%%      
total_num=0;
[size_x size_y]=size(output);
for i=1:size_x
    for j=1:size_y
        if output(i,j)>10^(-4)
            total_num=total_num+1;
        end
    end
end
clear size_x size_y i j
%%
align_data=zeros(0.0);
counting=1;
j=1;

for i=1:total_num
    if mod(counting,10001)==0
        j=j+1;
        counting=1;
        index=find(temp_output==max(temp_output(:)));
        [x y]=size(index);
         for k=1:x
           align_data(counting,j)=temp_output(index(k));
           temp_output(index(k))=0;
           counting=counting+1;
         end
    else 
        index=find(temp_output==max(temp_output(:)));
        [x y]=size(index);
         for k=1:x
           align_data(counting,j)=temp_output(index(k));
           temp_output(index(k))=0;
           counting=counting+1;
         end
    end
end

clear x y i j k counting temp_output index