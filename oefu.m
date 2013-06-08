function [g,estavel]=oefu(u,y,ModelParam,FF);
%

    OPT.dsp=0;
    estavel=[];
    mmm=[];
    mmm.type='oe';
    mmm.nB=ModelParam(2)-1;
    mmm.nA=ModelParam(1);
    mmm.delay=ModelParam(3);

    mi=mmm;

for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(3,FF(i));
    yf=filter(fb,fa,y);
    uf=filter(fb,fa,u);  

    %Estimate with filtered data
    zzz.u=uf';
    zzz.y=yf';
    g=est(zzz,mmm,OPT);

    %Test if the model is stable
    if max(abs(roots(g.A)))>1
        mmm=mi;
        estavel(i)=1;
    else
        mmm=g;
        estavel(i)=0;
    end    
    
end


    %Estimate with real data
    zzz.u=u';
    zzz.y=y';
    g=est(zzz,mmm,OPT);
