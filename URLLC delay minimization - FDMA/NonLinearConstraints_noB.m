function f = NonLinearConstraints_noB(L,g,No,C,B,Eu,Em,Tm,U,M,K,scaleF)
    function [c,ceq] = mycon(x)
        Tu = x(1);
        Fu = scaleF.*x(1+1 : 1+U);
        pu = x(1+U+1 : 1+2*U);
        Fm = scaleF.*x(1+2*U+1 : 1+2*U+M);
        pm = x(1+2*U+M+1 : end);
                
        latency = zeros(K,1);
        energy = zeros(K,1);
        for i = 1:U
            latency(i)= L/(B*log2(1+(g(i)*pu(i))/(B*No))) + L*C/Fu(i) - Tu;
            latency(i+U)= L/(B*log2(1+(g(i+U)*pm(i))/(B*No))) + L*C/Fm(i) - Tm;
            energy(i) = L*pu(i)- Eu*(B*log2(1+(g(i)*pu(i))/(B*No)));
            energy(i+U) = L*pm(i) - Em*(B*log2(1+(g(i+U)*pm(i))/(B*No)));
        end
        c = [latency; energy];
        ceq = [];
       
    end
f = @mycon;
end