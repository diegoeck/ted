function [G]=simptf(G);
%This function simplify a transfer function G
%
%Usage: [G]=simptf(G)
%
%Diego Eckhard - 25/06/2014
%UFRGS Identification Toolbox
t=size(G);
k1=t(1);
k2=t(2);
for ii=1:k1
for jj=1:k2

    
[z,p,k,t]=zpkdata(G(ii,jj),'v');


a=0;
while a==0
    a=1;

    for i=1:length(p)
        
        for j=1:length(z)

            if round(1000*p(i))==round(1000*z(j))
                
                
                polo=p(i);
                zero=z(j);
                
                p(i)=[];
                z(j)=[];
                
             
                
                if imag(polo)>10^-4
                    p(i)=[];
                    z(j)=[];
%                     for iii=i:length(p)
%                        if round(1000*polo')==round(1000*p(iii))
%                         p(iii)=[];
%                         
%                        end
%                       break
%                     end
%                     for jjj=j:length(z)
%                        if round(1000*zero')==round(1000*z(jjj))
%                         z(jjj)=[];
%                        end
%                        break
%                     end
%                 
                end
                
                a=0;
                break
            end
        end
        if a==0
            break
        end

    end
end

for i=1:length(p)
if abs(imag(p(i)))<10^-4
    p(i)=real(p(i));
end
end

for i=1:length(z)
if abs(imag(z(i)))<10^-4
    z(i)=real(z(i));
end
end

if abs(k)<10^-4
    k=0
end

G(ii,jj)=zpk(z,p,k,t);




end
end