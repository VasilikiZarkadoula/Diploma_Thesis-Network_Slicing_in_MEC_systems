function f = NonLinearConstraints_noB_scaling(R,g,No,C,B,E,E2,T2,M,N,K,scaleF)
    function [c,ceq] = mycon(x)
        time_threshold = x(1);
        Fm = scaleF.*x(1+1 : 1+M);
        pm = x(1+M+1 : 1+2*M);
        Fn = scaleF.*x(1+2*M+1 : 1+2*M+N);
        pn = x(1+2*M+N+1 : end);
                
        latency = zeros(K,1);
        energy = zeros(K,1);
        for i = 1:M
            latency(i)= R/(B*log2(1+(g(i)*pm(i))/(B*No))) + R*C/Fm(i) - time_threshold;
            latency(i+M)= R/(B*log2(1+(g(i+M)*pn(i))/(B*No))) + R*C/Fn(i) - T2;
            energy(i) = R*pm(i)- E*(B*log2(1+(g(i)*pm(i))/(B*No)));
            energy(i+M) = R*pn(i) - E2*(B*log2(1+(g(i+M)*pn(i))/(B*No)));
        end
        c = [latency; energy];
        ceq = [];
       
    end
f = @mycon;
end