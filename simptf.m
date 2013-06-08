function [G]=simtf(G);
%This function simplify a transfer function G
%
%Usage: [G]=simtf(G)
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox
[z,p,k,t]=zpkdata(G,'v');

a=0;
while a==0
    a=1;

    for i=1:length(p)
        
        for j=1:length(z)

            if round(1000*p(i))==round(1000*z(j))
                p(i)=[];
                z(j)=[];
                a=0;
                break
            end
        end
        if a==0
            break
        end

    end
end
G=zpk(z,p,k,t);