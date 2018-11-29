%% from raw data
% l_aws_data = zeros(23,1);
% 
% data = importdata('1995.csv',',',1);
% aws_data = data.data;
% aws_data(~isfinite(aws_data))=0;
% l_aws_data(1,1) = length(aws_data);
% clear data
% 
% data = importdata('1996.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(2,1) = length(data);
% clear data
% 
% 
% data = importdata('1997.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(3,1) = length(data);
% clear data
% data = importdata('1998.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(4,1) = length(data);
% clear data
% data = importdata('1999.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(5,1) = length(data);
% clear data
% data = importdata('2000.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(6,1) = length(data);
% clear data
% data = importdata('2001.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(7,1) = length(data);
% clear data
% data = importdata('2002.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(8,1) = length(data);
% clear data
% data = importdata('2003.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(9,1) = length(data);
% clear data
% data = importdata('2004.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(10,1) = length(data);
% clear data
% data = importdata('2005.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(11,1) = length(data);
% clear data
% data = importdata('2006.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(12,1) = length(data);
% clear data
% data = importdata('2007.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(13,1) = length(data);
% clear data
% data = importdata('2008.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(14,1) = length(data);
% clear data
% data = importdata('2009.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(15,1) = length(data);
% clear data
% data = importdata('2010.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(16,1) = length(data);
% clear data
% data = importdata('2011.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(17,1) = length(data);
% clear data
% data = importdata('2012.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(18,1) = length(data);
% clear data
% data = importdata('2013.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(19,1) = length(data);
% clear data
% data = importdata('2014.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(20,1) = length(data);
% clear data
% data = importdata('2015.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(21,1) = length(data);
% clear data
% data = importdata('2016.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(22,1) = length(data);
% clear data
% data = importdata('2017.csv',',',1);
% data = data.data;
% data(~isfinite(data))=0;
% aws_data = [aws_data;data];
% l_aws_data(23,1) = length(data);
% clear data
% 


%%
n = length(aws_data);
for i = 1:n
    if mod(aws_data(i,2),1/24)<10/1440
        aws_data(i,2) = aws_data(i,2) - mod(aws_data(i,2),1/24);
    else
        aws_data(i,2) = aws_data(i,2) + (1/24 - mod(aws_data(i,2),1/24));
    end
end
clear i
result_t = unique(aws_data(:,2)); 
result_ws = result_t;
result_wd = result_t;
result_p = result_t;

result_t(1,21)=0;
result_ws(1,21)=0;
result_wd(1,21)=0;
result_p(1,21)=0;

% for i = 1:n
%     if aws_data(i,1)==400
%         result_t(find(aws_data(i,2)==result_t(:,1)),2) = aws_data(i,3);
%         result_wd(find(aws_data(i,2)==result_t(:,1)),2) = aws_data(i,4);
%         result_ws(find(aws_data(i,2)==result_t(:,1)),2) = aws_data(i,5);
%         result_p(find(aws_data(i,2)==result_t(:,1)),2) = aws_data(i,6);
%     end
% end
% clear i
for i = 1:n
    if aws_data(i,1)==401
        result_t(find(aws_data(i,2)==result_t(:,1)),3) = aws_data(i,3);
        result_wd(find(aws_data(i,2)==result_t(:,1)),3) = aws_data(i,4);
        result_ws(find(aws_data(i,2)==result_t(:,1)),3) = aws_data(i,5);
        result_p(find(aws_data(i,2)==result_t(:,1)),3) = aws_data(i,6);
    end
end
clear i




% 
% 
% 
% 
% for i = 1:n
%     if aws_data(i,1)==410
%         aws_data(i,:)=[];
%         i=i-1;
%         n=n-1;
%     elseif aws_data(i,1)==417
%         aws_data(i,:)=[];
%         i=i-1;
%         n=n-1;
%     elseif aws_data(i,1)==421
%         aws_data(i,:)=[];
%         i=i-1;
%         n=n-1;
%     elseif aws_data(i,1)==422
%         aws_data(i,:)=[];
%         i=i-1;
%         n=n-1;
%     elseif aws_data(i,1)==423
%         aws_data(i,:)=[];
%         i=i-1;
%         n=n-1;
%     elseif aws_data(i,1)==424
%         aws_data(i,:)=[];
%         i=i-1;
%         n=n-1;
%     elseif aws_data(i,1)==425
%         aws_data(i,:)=[];
%         i=i-1;
%         n=n-1;
%     elseif aws_data(i,1)==889
%         aws_data(i,:)=[];
%         i=i-1;
%         n=n-1;
%     end
% end
% clear i
% result2 = unique(aws_data(:,1));
% 




