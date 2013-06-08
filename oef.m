function [Model,estavel]=oef(u,y,order,FF);
%

estavel=[];


     [fb,fa] = butter(3,FF(1));
     yf=filter(fb,fa,y);
     uf=filter(fb,fa,u); 

dat=iddata(yf,uf,1);
Model = arx(dat,order);
[mA,mB,mC,mD,mF] = polydata(Model);
Model = idpoly(1,mB,1,1,mA,1,1);


for i=1:length(FF)
    i;
    
     [fb,fa] = butter(3,FF(i));
     yf=filter(fb,fa,y);
     uf=filter(fb,fa,u);  

dat=iddata(yf,uf,1);
Model = oe(dat,Model);

[mA,mB,mC,mD,mF] = polydata(Model);
 

if max(abs(roots(mF)))>1
dat=iddata(yf,uf,1);
Model = arx(dat,order);
[mA,mB,mC,mD,mF] = polydata(Model);
Model = idpoly(1,mB,1,1,mA,1,1);
    if max(abs(roots(mA)))>1
        Model=order;
    end

    estavel(i)=1;
else
    estavel(i)=0;
    Modele=Model;
end    
    
    
end

max(abs(roots(mF)));

dat=iddata(y,u,1);
Model = oe(dat,Model);
