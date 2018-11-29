n_data = [1508568; 1687032; 1830840; 1940376; 2120016; 2190336; 2302896; 2444939; 2496648; 2543111; 2635183; 2608779; 2766000];
n_data2 = zeros(13,2);
% 	% 1: 111152, 2: 823691 
n_stn1 = 823651;
n_stn2 = 823652;
% 
% for i = 1:n_data(1,1)
% 		if x2001(i,1) == n_stn1
%            n_data2(1,1) =  n_data2(1,1)+1;
%         elseif x2001(i,1) == n_stn2
%            n_data2(1,2) =  n_data2(1,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(2,1)
% 		if x2002(i,1) == n_stn1
%            n_data2(2,1) =  n_data2(2,1)+1;
%         elseif x2002(i,1) == n_stn2
%            n_data2(2,2) =  n_data2(2,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(3,1)
% 		if x2003(i,1) == n_stn1
%            n_data2(3,1) =  n_data2(3,1)+1;
%         elseif x2003(i,1) == n_stn2
%            n_data2(3,2) =  n_data2(3,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(4,1)
% 		if x2004(i,1) == n_stn1
%            n_data2(4,1) =  n_data2(4,1)+1;
%         elseif x2004(i,1) == n_stn2
%            n_data2(4,2) =  n_data2(4,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(5,1)
% 		if x2005(i,1) == n_stn1
%            n_data2(5,1) =  n_data2(5,1)+1;
%         elseif x2005(i,1) == n_stn2
%            n_data2(5,2) =  n_data2(5,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(6,1)
% 		if x2006(i,1) == n_stn1
%            n_data2(6,1) =  n_data2(6,1)+1;
%         elseif x2006(i,1) == n_stn2
%            n_data2(6,2) =  n_data2(6,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(7,1)
% 		if x2007(i,1) == n_stn1
%            n_data2(7,1) =  n_data2(7,1)+1;
%         elseif x2007(i,1) == n_stn2
%            n_data2(7,2) =  n_data2(7,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(8,1)
% 		if x2008(i,1) == n_stn1
%            n_data2(8,1) =  n_data2(8,1)+1;
%         elseif x2008(i,1) == n_stn2
%            n_data2(8,2) =  n_data2(8,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(9,1)
% 		if x2009(i,1) == n_stn1
%            n_data2(9,1) =  n_data2(9,1)+1;
%         elseif x2009(i,1) == n_stn2
%            n_data2(9,2) =  n_data2(9,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(10,1)
% 		if x2010(i,1) == n_stn1
%            n_data2(10,1) =  n_data2(10,1)+1;
%         elseif x2010(i,1) == n_stn2
%            n_data2(10,2) =  n_data2(10,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(11,1)
% 		if x2011(i,1) == n_stn1
%            n_data2(11,1) =  n_data2(11,1)+1;
%         elseif x2011(i,1) == n_stn2
%            n_data2(11,2) =  n_data2(11,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(12,1)
% 		if x2012(i,1) == n_stn1
%            n_data2(12,1) =  n_data2(12,1)+1;
%         elseif x2012(i,1) == n_stn2
%            n_data2(12,2) =  n_data2(12,2)+1;
%         end
% end
% clear i
% 
% for i = 1:n_data(13,1)
% 		if x2013(i,1) == n_stn1
%            n_data2(13,1) =  n_data2(13,1)+1;
%         elseif x2013(i,1) == n_stn2
%            n_data2(13,2) =  n_data2(13,2)+1;
%         end
% end
% clear i
% 

n=0;
for i = 1:13
    n = n + n_data2(i,1);
end
clear i
stn1 = zeros(n,7);

n=0;
for i = 1:13
    n = n + n_data2(i,2);
end
clear i
stn2 = zeros(n,7);
clear n

n1=0;
n2=0;
for i = 1:n_data(1,1)
		if x2001(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2001(i,j);
            end
            clear j
        elseif x2001(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2001(i,j);
            end
            clear j            
        end
end
clear i
for i = 1:n_data(2,1)
		if x2002(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2002(i,j);
            end
            clear j           
        elseif x2002(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2002(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(3,1)
		if x2003(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2003(i,j);
            end
            clear j           
        elseif x2003(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2003(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(4,1)
		if x2004(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2004(i,j);
            end
            clear j           
        elseif x2004(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2004(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(5,1)
		if x2005(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2005(i,j);
            end
            clear j           
        elseif x2005(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2005(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(6,1)
		if x2006(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2006(i,j);
            end
            clear j           
        elseif x2006(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2006(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(7,1)
		if x2007(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2007(i,j);
            end
            clear j           
        elseif x2007(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2007(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(8,1)
		if x2008(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2008(i,j);
            end
            clear j           
        elseif x2008(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2008(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(9,1)
		if x2009(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2009(i,j);
            end
            clear j           
        elseif x2009(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2009(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(10,1)
		if x2010(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2010(i,j);
            end
            clear j           
        elseif x2010(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2010(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(11,1)
		if x2011(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2011(i,j);
            end
            clear j                      
        elseif x2011(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2011(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(12,1)
		if x2012(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2012(i,j);
            end
            clear j                      
        elseif x2012(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2012(i,j);
            end
            clear j                       
        end
end
clear i

for i = 1:n_data(13,1)
		if x2013(i,1) == n_stn1
            n1 = n1+1;
            for j = 1:7
                stn1(n1,j) = x2013(i,j);
            end
            clear j                      
        elseif x2013(i,1) == n_stn2
            n2 = n2 +1;
            for j = 1:7
                stn2(n2,j) = x2013(i,j);
            end
            clear j                       
        end
end
clear i
