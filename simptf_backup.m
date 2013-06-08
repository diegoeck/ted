function [G]=simtf(G);
%This function simplify a transfer function G
%
%Usage: [G]=simtf(G)
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox
[z,p,kk,t]=zpkdata(G,'v');
zz=z;
pp=p;

for i=1:length(p)

k=0;
while k<length(pp)
    k=k+1;
    l=0;
    while l<length(zz)
        l=l+1;
         if round(1000*pp(k))==round(1000*zz(l))
             pp(k)=[];
             zz(l)=[];
         end
    end
end

end
G=zpk(zz,pp,kk,t);