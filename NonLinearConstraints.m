function f = NonLinearConstraints(R,h,No,C,E,E2,T2,M,N)
    function [c,ceq] = mycon(x)
        % low latency users constraints
        m_latency = zeros(M,1);
        m_latency_offload = zeros(M,1);
        m_energy = zeros(M,1);
        m_energy_offload = zeros(M,1);
        for i = 0:M-1
            % constraint 6
            m_latency(i+1) = (R(i+1))/(x(2+i)*log(1+(x(2+2*M+i)*(h)^2)/(x(2+i)*No))) + (C(i+1)*R(i+1)/x(2+M+i)) - x(1);
            m_latency_offload(i+1) = -((R(i+1)) / (x(2+i)*log(1+(x(2+2*M+i)*(h)^2)/(x(2+i)*No))) + (C(i+1)*R(i+1))/x(2+M+i));
            % constraint 8
            m_energy(i+1) = (R(i+1)*x(2+2*M+i))/(x(2+i)*log(1+(x(2+2*M+i)*(h)^2)/(x(2+i)*No))) - E;
            m_energy_offload(i+1) = -(R(i+1)*x(2+2*M+i))/(x(2+i)*log(1+(x(2+2*M+i)*(h)^2)/(x(2+i)*No)));
        end
        
        % low energy consumption users constraints
        n_latency = zeros(N,1);
        n_latency_offload = zeros(N,1);
        n_energy = zeros(N,1);
        n_energy_offload = zeros(N,1);
        for j = 0:N-1
            % constraint 7
            n_latency(j+1) = (R(M+j+1))/(x(2+3*M+j)*log(1+(x(2+3*M+2*N+j)*(h)^2)/(x(2+3*M+j)*No))) + (C(M+j+1)*R(M+j+1))/x(2+3*M+N+j) - T2;
            n_latency_offload(j+1) = -((R(M+j+1))/(x(2+3*M+j)*log(1+(x(2+3*M+2*N+j)*(h)^2)/(x(2+3*M+j)*No))) + (C(M+j+1)*R(M+j+1))/x(2+3*M+N+j));
            % constraint 9
            n_energy(j+1) = (R(M+j+1)*x(2+3*M+2*N+j))/(x(2+3*M+j)*log(1+(x(2+3*M+2*N+j)*(h)^2)/(x(2+3*M+j)*No))) - E2;
            n_energy_offload(j+1) = -((R(M+j+1)*x(2+3*M+2*N+j))/(x(2+3*M+j)*log(1+(x(2+3*M+2*N+j)*(h)^2)/(x(2+3*M+j)*No))));
        end
        
        c = [m_latency; m_latency_offload; m_energy; m_energy_offload;...
            n_latency; n_latency_offload; n_energy; n_energy_offload];
        ceq = [];
    end
f = @mycon;
end