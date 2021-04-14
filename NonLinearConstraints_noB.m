function f = NonLinearConstraints_noB(R,g,No,C,B,E,E2,T2,M,N)
    function [c,ceq] = mycon(x)
        % low latency users constraints
        m_latency = zeros(M,1);
        m_latency_offload = zeros(M,1);
        m_energy = zeros(M,1);
        m_energy_offload = zeros(M,1);
        for i = 1:M
            % constraint 6
            m_latency(i) = (R(i))/(B(i)*log2(1+(x(1+M+i)*g(i))/(B(i)*No))) + (C(i)*R(i))/x(1+i) - x(1);
            m_latency_offload(i) = -((R(i))/(B(i)*log2(1+(x(1+M+i)*g(i))/(B(i)*No))) + (C(i)*R(i))/x(1+i));
            % constraint 8
            m_energy(i) = (R(i)*x(1+M+i))/(B(i)*log2(1+(x(1+M+i)*g(i))/(B(i)*No))) - E;
            m_energy_offload(i) = -((R(i)*x(1+M+i))/(B(i)*log2(1+(x(1+M+i)*g(i))/(B(i)*No))));
        end
        
        % low energy consumption users constraints
        n_latency = zeros(N,1);
        n_latency_offload = zeros(N,1);
        n_energy = zeros(N,1);
        n_energy_offload = zeros(N,1);
        for j = 1:N
            % constraint 7
            n_latency(j) = (R(M+j))/(B(M+j)*log2(1+(x(1+2*M+N+j)*g(M+j))/(B(M+j)*No))) + (C(M+j)*R(M+j))/x(1+2*M+j) - T2;
            n_latency_offload(j) = -((R(M+j))/(B(M+j)*log2(1+(x(1+2*M+N+j)*g(M+j))/(B(M+j)*No))) + (C(M+j)*R(M+j))/x(1+2*M+j));
            % constraint 9
            n_energy(j) = (R(M+j)*x(1+2*M+N+j))/(B(M+j)*log2(1+(x(1+2*M+N+j)*g(M+j))/(B(M+j)*No))) - E2;
            n_energy_offload(j) = -((R(M+j)*x(1+2*M+N+j))/(B(M+j)*log2(1+(x(1+2*M+N+j)*g(M+j))/(B(M+j)*No))));
        end
        
        c = [m_latency; m_latency_offload; m_energy; m_energy_offload;...
            n_latency; n_latency_offload; n_energy; n_energy_offload];
        ceq = [];
    end
f = @mycon;
end