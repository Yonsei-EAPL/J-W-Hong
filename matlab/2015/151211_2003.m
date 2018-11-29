% G = map2mat(ones(size(H,2),size(H,3)),H);
% 
% [E,pc,expvar] = caleof(G,N,method);  
% eof = mat2map(ones(size(H,2),size(H,3)),E);
% 
% plot(seoul(:,1),seoul(:,2),'r')

%% data extract (only near by Seoul)
% raw2003 = aws2003;
% temp = 0;
% length=3711259;
% for i = 1:length
%     if aws2003(i,3)<37
%         temp = temp+1;
%     elseif aws2003(i,4)>127.5
%         temp = temp+1;
%     end
% end
% clear i 
% aws2003 = zeros(length-temp,12);
% temp=0;
% for i = 1:length
%     if (raw2003(i,3)>37)&&(raw2003(i,4)<127.5)
%         temp = temp+1;
%         aws2003(temp,:) = raw2003(i,:);
%     end
% end
% clear i temp


%% site information
length = 709228;
temp = 0;
infs = 0;
for i = 1:length
    if i ==1
        infs = aws2003(i,2);
        temp = temp+1;
    elseif infs~=aws2003(i,2)
        infs = aws2003(i,2);
        temp = temp+1;
    end
end
clear i infs
siteinf = zeros(83,3);
temp = 0;
infs = 0;
for i = 1:length
    if i ==1
        infs = aws2003(i,2);
        temp = temp+1;
        siteinf(temp,1) = infs;
        siteinf(temp,2) = aws2003(i,3);
        siteinf(temp,3) = aws2003(i,4);
    elseif infs~=aws2003(i,2)
        infs = aws2003(i,2);
        temp = temp+1;
        siteinf(temp,1) = infs;
        siteinf(temp,2) = aws2003(i,3);
        siteinf(temp,3) = aws2003(i,4);        
    end
end
clear i infs temp


%%


