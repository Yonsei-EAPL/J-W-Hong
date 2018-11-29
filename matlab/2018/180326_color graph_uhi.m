x=0
y=0
a=0
c=0
s1=0
s2=0



s2=a;
a=a*100;
s1=num2str(s1);
s2=round(s2,1);
s2=num2str(s2);

figure()
scatter(x,y,a,c,'filled')
hold on

text(x,y,s1)
%text(x,y-10,s2)
text(x,y-10,s2)