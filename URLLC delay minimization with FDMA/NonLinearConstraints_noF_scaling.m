function f = NonLinearConstraints_noF_scaling(R,g,No,C,F,E,E2,T2,M,N,K,scaleB)
    function [c,ceq] = mycon(x)
        time_threshold = x(1);
        Bm = scaleB.*x(1+1 : 1+M);
        pm = x(1+M+1 : 1+2*M);
        Bn = scaleB.*x(1+2*M+1 : 1+2*M+N);
        pn = x(1+2*M+N+1 : end);
                
        latency = zeros(K,1);
        energy = zeros(K,1);
        for i = 1:M
            latency(i)= R/(Bm(i)*log2(1+(g(i)*pm(i))/(Bm(i)*No))) + R*C/F - time_threshold;
            latency(i+M)= R/(Bn(i)*log2(1+(g(i+M)*pn(i))/(Bn(i)*No))) + R*C/F - T2;
            energy(i) = R*pm(i)- E*(Bm(i)*log2(1+(g(i)*pm(i))/(Bm(i)*No)));
            energy(i+M) = R*pn(i) - E2*(Bn(i)*log2(1+(g(i+M)*pn(i))/(Bn(i)*No)));
        end
        c = [latency; energy];
        ceq = [];
       
    end
f = @mycon;
end