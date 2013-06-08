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

    zzz.u=u';
    zzz.y=y';
%     ma.type='arx';
%     ma=est(zzz,ma,OPT);
%     ma.type='oe';
%     if max(abs(roots(ma.A)))>1
%         ma=mi;
%     end
    
    

for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(5,FF(i));
    yf=filter(fb,fa,y);
    uf=filter(fb,fa,u);  

    %Estimate with filtered data
    zzz.u=uf';
    zzz.y=yf';
    g=est(zzz,mmm,OPT);

    %Test if the model is stable
    if max(abs(roots(g.A)))>1

            m2=mi;
            m2.type='arx';
            m2=est(zzz,m2,OPT);
            m2.type='oe';
            mmm=m2;
            estavel(i)=1;
            
            if max(abs(roots(m2.A)))>1
                mmm=mi;
                estavel(i)=2;
            end
    else
        mmm=g;
        estavel(i)=0;
    end    
    
end


    %Estimate with real data
    zzz.u=u';
    zzz.y=y';
    g=est(zzz,mmm,OPT);
