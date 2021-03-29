function f = constraints6to9(r,h,no,cy,e,e2,t2,m,n)
    function [c,ceq] = mycon(x)
        % low latency users constraints
        m_latency = zeros(m,1);
        m_latency_offload = zeros(m,1);
        m_energy = zeros(m,1);
        m_energy_offload = zeros(m,1);
        for i = 0:m-1
            % constraint 6
            m_latency(i+1) = r/(x(2+i)*log(1+(x(2+2*m+i)*(h)^2)/(x(2+i)*no))) + cy*r/x(2+m+i) - x(1);
            m_latency_offload(i+1) = -(r/(x(2+i)*log(1+(x(2+2*m+i)*(h)^2)/(x(2+i)*no))) + cy*r/x(2+m+i));
            % constraint 8
            m_energy(i+1) = (r*x(2+2*m+i))/(x(2+i)*log(1+(x(2+2*m+i)*(h)^2)/(x(2+i)*no))) - e;
            m_energy_offload(i+1) = -(r*x(2+2*m+i))/(x(2+i)*log(1+(x(2+2*m+i)*(h)^2)/(x(2+i)*no)));
        end
        
        % low energy consumption users constraints
        n_latency = zeros(n,1);
        n_latency_offload = zeros(n,1);
        n_energy = zeros(n,1);
        n_energy_offload = zeros(n,1);
        for j = 0:n-1
            % constraint 7
            n_latency(j+1) = r/(x(2+3*m+j)*log(1+(x(2+3*m+2*n+j)*(h)^2)/(x(2+3*m+j)*no))) + cy*r/x(2+3*m+n+j) - t2;
            n_latency_offload(j+1) = -(r/(x(2+3*m+j)*log(1+(x(2+3*m+2*n+j)*(h)^2)/(x(2+3*m+j)*no))) + cy*r/x(2+3*m+n+j));
            % constraint 9
            n_energy(j+1) = (r*x(2+3*m+2*n+j))/(x(2+3*m+j)*log(1+(x(2+3*m+2*n+j)*(h)^2)/(x(2+3*m+j)*no))) - e2;
            n_energy_offload(j+1) = -((r*x(2+3*m+2*n+j))/(x(2+3*m+j)*log(1+(x(2+3*m+2*n+j)*(h)^2)/(x(2+3*m+j)*no))));
        end
        
        c = [m_latency; m_latency_offload; m_energy; m_energy_offload;...
            n_latency; n_latency_offload; n_energy; n_energy_offload];
        ceq = [];
    end
f = @mycon;
end