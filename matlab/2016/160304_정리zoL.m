zol0=0;
zol01=0;
zol03=0;
zol05=0;
zol10=0;
zol20=0;
zol00=0;

l=38909;
for i=1:l
    if ((zol(i,1)>-2.0)&&(zol(i,1)<-0.1))
        zol00= zol00+1;
        zol(i,2)=-2;
    elseif ((zol(i,1)>-0.01)&&(zol(i,1)<0.01))
        zol0= zol0+1;
        zol(i,2)=0;
    elseif ((zol(i,1)>0.09)&&(zol(i,1)<0.11))
        zol01= zol01+1;
        zol(i,2)=0.1;
    elseif ((zol(i,1)>0.29)&&(zol(i,1)<0.31))
        zol03= zol03+1;
        zol(i,2)=0.3;
    elseif ((zol(i,1)>0.49)&&(zol(i,1)<0.51))
        zol05= zol05+1;
        zol(i,2)=0.5;
    elseif ((zol(i,1)>0.9)&&(zol(i,1)<1.1))
        zol10= zol10+1;
        zol(i,2)=1.0;
    elseif ((zol(i,1)>1.9)&&(zol(i,1)<2.1))
        zol20= zol20+1;
        zol(i,2)=2.0;
    else
        zol(i,2)=-999;
    end
end
clear i