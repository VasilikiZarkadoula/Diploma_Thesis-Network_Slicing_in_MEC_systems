function f = NonLinearConstraints_noF(L,g,No,C,F,Eu,Em,Tm,U,M,K,scaleB)
    function [c,ceq] = mycon(x)
        Tu = x(1);
        Bu = scaleB.*x(1+1 : 1+U);
        pu = x(1+U+1 : 1+2*U);
        Bm = scaleB.*x(1+2*U+1 : 1+2*U+M);
        pm = x(1+2*U+M+1 : end);
                
        latency = zeros(K,1);
        energy = zeros(K,1);
        for i = 1:U
            latency(i)= L/(Bu(i)*log2(1+(g(i)*pu(i))/(Bu(i)*No))) + L*C/F - Tu;
            latency(i+U)= L/(Bm(i)*log2(1+(g(i+U)*pm(i))/(Bm(i)*No))) + L*C/F - Tm;
            energy(i) = L*pu(i)- Eu*(Bu(i)*log2(1+(g(i)*pu(i))/(Bu(i)*No)));
            energy(i+U) = L*pm(i) - Em*(Bm(i)*log2(1+(g(i+U)*pm(i))/(Bm(i)*No)));
        end
        c = [latency; energy];
        ceq = [];
       
    end
f = @mycon;
end