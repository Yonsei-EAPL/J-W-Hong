% 108 site

for i = 1:length(x)
   if (x(i,2)<0)||(x(i,2)>360)
       x(i,3) = 888.8;
   elseif x(i,4)<2013010000
       if x(i,2)<45
           x(i,3) = 18;
       elseif x(i,2)<90
           x(i,3) = 138;
       elseif x(i,2)<135
           x(i,3) = 138;
       elseif x(i,2)<180
           x(i,3) = 69;
       elseif x(i,2)<225
           x(i,3) = 18;
       elseif x(i,2)<270
           x(i,3) = 18;
       elseif x(i,2)<315
           x(i,3) = 18;
       else
           x(i,3) = 18;
       end
   elseif x(i,4)<2017010000
       if x(i,2)<45
           x(i,3) = 18;
       elseif x(i,2)<90
           x(i,3) = 138;
       elseif x(i,2)<135
           x(i,3) = 138;
       elseif x(i,2)<180
           x(i,3) = 69;
       elseif x(i,2)<225
           x(i,3) = 223;
       elseif x(i,2)<270
           x(i,3) = 223;
       elseif x(i,2)<315
           x(i,3) = 223;
       else
           x(i,3) = 18;
       end
   else
       if x(i,2)<45
           x(i,3) = 18;
       elseif x(i,2)<90
           x(i,3) = 138;
       elseif x(i,2)<135
           x(i,3) = 138;
       elseif x(i,2)<180
           x(i,3) = 69;
       elseif x(i,2)<225
           x(i,3) = 52;
       elseif x(i,2)<270
           x(i,3) = 52;
       elseif x(i,2)<315
           x(i,3) = 52;
       else
           x(i,3) = 18;
       end
    end    
end
clear i
