function f = NonLinearConstraints(L,g,No,C,Eu,Em,Tm,U,M,K,scaleB,scaleF)
    function [c,ceq] = mycon(x)
        Tu = x(1);
        Bu = scaleB.*x(1+1 : 1+U);
        Fu = scaleF.*x(1+U+1 : 1+2*U);
        pu = x(1+2*U+1 : 1+3*U);
        Bm = scaleB.*x(1+3*U+1 : 1+3*U+M);
        Fm = scaleF.*x(1+3*U+M+1 : 1+3*U+2*M);
        pm = x(1+3*U+2*M+1 : end);
                
        latency = zeros(K,1);
        energy = zeros(K,1);
        for i = 1:U
            latency(i)= L/(Bu(i)*log2(1+(g(i)*pu(i))/(Bu(i)*No))) + L*C/Fu(i) - Tu;
            latency(i+U)= L/(Bm(i)*log2(1+(g(i+U)*pm(i))/(Bm(i)*No))) + L*C/Fm(i) - Tm;
            energy(i) = L*pu(i)- Eu*(Bu(i)*log2(1+(g(i)*pu(i))/(Bu(i)*No)));
            energy(i+U) = L*pm(i) - Em*(Bm(i)*log2(1+(g(i+U)*pm(i))/(Bm(i)*No)));
        end
        c = [latency; energy];
        ceq = [];
       
    end
f = @mycon;
end