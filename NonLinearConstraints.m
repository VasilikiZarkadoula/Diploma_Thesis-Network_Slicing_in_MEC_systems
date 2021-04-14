function f = NonLinearConstraints(R,g,No,C,E,E2,T2,M,N)
    function [c,ceq] = mycon(x)
        % low latency users constraints
        m_latency = zeros(M,1);
        m_latency_offload = zeros(M,1);
        m_energy = zeros(M,1);
        m_energy_offload = zeros(M,1);
        for i = 1:M
            % constraint 6
            m_latency(i) = (R(i))/(x(1+i)*log2(1+(x(1+2*M+i)*g(i))/(x(1+i)*No))) + (C(i)*R(i)/x(1+M+i)) - x(1);
            m_latency_offload(i) = -((R(i)) / (x(1+i)*log2(1+(x(1+2*M+i)*g(i))/(x(1+i)*No))) + (C(i)*R(i))/x(1+M+i));
            % constraint 8
            m_energy(i) = (R(i)*x(1+2*M+i))/(x(1+i)*log2(1+(x(1+2*M+i)*g(i))/(x(1+i)*No))) - E;
            m_energy_offload(i) = -((R(i)*x(1+2*M+i))/(x(1+i)*log2(1+(x(1+2*M+i)*g(i))/(x(1+i)*No))));
        end
        
        % low energy consumption users constraints
        n_latency = zeros(N,1);
        n_latency_offload = zeros(N,1);
        n_energy = zeros(N,1);
        n_energy_offload = zeros(N,1);
        for j = 1:N
            % constraint 7
            n_latency(j) = (R(M+j))/(x(1+3*M+j)*log2(1+(x(1+3*M+2*N+j)*g(M+j))/(x(1+3*M+j)*No))) + (C(M+j)*R(M+j))/x(1+3*M+N+j) - T2;
            n_latency_offload(j) = -((R(M+j))/(x(1+3*M+j)*log2(1+(x(1+3*M+2*N+j)*g(M+j))/(x(1+3*M+j)*No))) + (C(M+j)*R(M+j))/x(1+3*M+N+j));
            % constraint 9
            n_energy(j) = (R(M+j)*x(1+3*M+2*N+j))/(x(1+3*M+j)*log2(1+(x(1+3*M+2*N+j)*g(M+j))/(x(2+3*M+j)*No))) - E2;
            n_energy_offload(j) = -((R(M+j)*x(1+3*M+2*N+j))/(x(1+3*M+j)*log2(1+(x(1+3*M+2*N+j)*g(M+j))/(x(1+3*M+j)*No))));
        end
        
        c = [m_latency; m_latency_offload; m_energy; m_energy_offload;...
            n_latency; n_latency_offload; n_energy; n_energy_offload];
        ceq = [];
    end
f = @mycon;
end