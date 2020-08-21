function [P]=curve(z,D)
    p=[0 z;0 0;D 0];
    n=3;
    n1=2;
    for    i=0:2
        sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  % for calculating (x!/(y!(x-y)!)) values 
    end
    l=[];
    UB=[];
    for u=0:0.2:1
        for d=1:3
           UB(d)=sigma(d)*((1-u)^(3-d))*(u^(d-1));
        end
        l=cat(1,l,UB);                                      %catenation 
    end
    P=l*p;
end
