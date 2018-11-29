%% split 2004
wait = waitbar(0,'...split...');
data = sfc_tim_2004;
for i = 1:676369
    waitbar((i/(676369*3)),wait,sprintf('%f',i/(676369*3)*100))
    data(i,1) = strrep(data(i,1),'=','-999');
    data(i,1) = strrep(data(i,1),'|',',');
    data(i,1) = strrep(data(i,1),',,',',-999,');
end
for i = 1:676369
    waitbar(((i+676369)/(676369*3)),wait,sprintf('%f',(i+676369)/(676369*3)*100))
    data(i,1) = strrep(data(i,1),',,',',');
end
data = regexp(data,',','split');
aws2004 = zeros(676368,11);
for i = 2:676369
    waitbar(((i+676369*2)/(676369*3)),wait,sprintf('%f',(i+676369*2)/(676369*3)*100))
    for j = 1:11
        temp1 = data{i,1}(1,j);
        aws2004(i-1,j) = str2num(temp1{1,1});
    end
end
clear i j temp1 data
close(wait);

%% statistics
wait = waitbar(0,'...statistics...');
temp = zeros(1,2);
for i = 1:676368
    waitbar((i/(676368)),wait,sprintf('%f',i/(676368)*100))
    [size_temp1 size_temp2] = size(temp);
    if temp(1,1) == 0
        temp(1,1) = aws2004(i,2);
        temp(1,2) = 1;
    else
        temp2 = 0;
        for j = 1:size_temp1
            if temp(j,1) == aws2004(i,2)
                temp2 = 1;
                temp(j,2) = temp(j,2)+1;
            end
        end
        if temp2 ==0
            temp(size_temp1+1,1) = aws2004(i,2);
            temp(size_temp1+1,2) = 1;
        end
    end
end
result = temp;
clear i j temp2 size_temp1 size_temp2 temp
close(wait);
