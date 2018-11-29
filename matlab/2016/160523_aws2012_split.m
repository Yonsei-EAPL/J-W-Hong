%% 2016-05-23 
% for ½ÉÈ­¿¬±¸
% with 2012.mat

stn138 = zeros(8768,5); % Æ÷Ç×1
stn253 = zeros(8784,5); % ±èÇØ2
stn371 = zeros(8662,5); % ±âÈï±¸3
stn406 = zeros(8780,5); % µµºÀ4
stn411 = zeros(8641,5); % ¸¶Æ÷5

temp1 = 0;
temp2 = 0;
temp3 = 0;
temp4 = 0;
temp5 = 0;

for i = 1:5981372
    if aws2012(i,2)==138
        temp1 = temp1+1;
        stn138(temp1,1) = aws2012(i,1);
        stn138(temp1,2) = aws2012(i,3);
        stn138(temp1,3) = aws2012(i,6);
        stn138(temp1,4) = aws2012(i,9);
        stn138(temp1,5) = aws2012(i,10);
    elseif aws2012(i,2)==253
        temp2 = temp2+1;
        stn253(temp2,1) = aws2012(i,1);
        stn253(temp2,2) = aws2012(i,3);
        stn253(temp2,3) = aws2012(i,6);
        stn253(temp2,4) = aws2012(i,9);
        stn253(temp2,5) = aws2012(i,10);
    elseif aws2012(i,2)==371
        temp3 = temp3+1;
        stn371(temp3,1) = aws2012(i,1);
        stn371(temp3,2) = aws2012(i,3);
        stn371(temp3,3) = aws2012(i,6);
        stn371(temp3,4) = aws2012(i,9);
        stn371(temp3,5) = aws2012(i,10);
    elseif aws2012(i,2)==406
        temp4 = temp4+1;
        stn406(temp4,1) = aws2012(i,1);
        stn406(temp4,2) = aws2012(i,3);
        stn406(temp4,3) = aws2012(i,6);
        stn406(temp4,4) = aws2012(i,9);
        stn406(temp4,5) = aws2012(i,10);
    elseif aws2012(i,2)==411
        temp5 = temp5+1;
        stn411(temp5,1) = aws2012(i,1);
        stn411(temp5,2) = aws2012(i,3);
        stn411(temp5,3) = aws2012(i,6);
        stn411(temp5,4) = aws2012(i,9);
        stn411(temp5,5) = aws2012(i,10);
    end
end
clear i temp1 temp2 temp3 temp4 temp5




