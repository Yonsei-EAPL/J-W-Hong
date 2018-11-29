x_data=zeros(max_tau+1,1);
z=5;

for i=1:max_tau+1
    x_data(i,1)=i-1;
end

num_figure_tot=0;
num_figure=0;
for i=1:num_var-1
    for k=i:num_var

        num_plot=(num_var-i)/z;
        a1=mod((num_var-i),z);
        if a1~=0
            num_plot=num_plot-a1/z+1;
        end
        
        for j=1:num_plot
            if j<numplot
                num_figure=num_figure+1;
                figure(num_figure);
                for k=1:5
                    
                    
                    
                    if ((j-1)*5+k)<(num_var-i)
                        y_max=zeros(4,1);
                        y_max(1,1)=max(graph1(:,((i-1)*num_var+(i+((j-1)*5+k)))));
                        y_max(2,1)=max(graph2(:,((i-1)*num_var+(i+((j-1)*5+k)))));
                        y_max(3,1)=max(graph1(:,(((j-1)*5+k)-1)*num_var+i));
                        y_max(4,1)=max(graph2(:,(((j-1)*5+k)-1)*num_var+i));
                        y_max_result=max(y_max);
                        subplot(z,1,k)
                        plot(x_data(:,1),graph1(:,((i-1)*num_var+(i+((j-1)*5+k)))))
                        axis([0 36 0 y_max_result])
                        hold on
                        plot(x_data(:,1),graph2(:,((i-1)*num_var+(i+((j-1)*5+k)))))
                        hold on
                        plot(x_data(:,1),graph1(:,(((j-1)*5+k)-1)*num_var+i))
                        hold on
                        plot(x_data(:,1),graph2(:,(((j-1)*5+k)-1)*num_var+i))
                        hold on
                    end
                end
            else
                
            end
        end
    end
end

