function [u]=experiment_realization_cos_par(r,tsinal);
%This function computes optimal input signal
%
%Usage: [u]=experiment_realization_cos_par(C,length);
%
%               C: optimal spectrum [c_1 c_2 ... c_N]
%
%               length: length of the input signal
%
%Diego Eckhard - 20/09/2011
%UFRGS Identification Toolbox



    N=length(r);
    sinal=randn(1,tsinal);


    
    r0=r(1);
    r1n=r(2:N)/2;
    
    A=toeplitz([r0 r1n(1:(N-2)) ]);
    
    %A=A(:,1:4)
    
    %eig(A)
    
    %Ai=pinv(A,0.01);
    
    %size(Ai)
    
    Ai=inv(A);
    f=Ai*r1n'
    
    ga=1+f'*r1n'
    
    %eig(A)
    
    %b=inv(sqrt(f'*r1n'+1))
    

    %abs(roots([1;-f]'))
    
    ro=roots([1;-f])
    abs(ro)
    acos(real(ro))
    
    %b=1
    u=filter([1],[1;-f]',sinal); 
    u=u';

    


