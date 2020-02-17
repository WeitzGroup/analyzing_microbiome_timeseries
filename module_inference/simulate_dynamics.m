function [t,H,V] = simulate_dynamics(params,x0,tfinal)
    [K,r,m,M,phi,beta] = params{:};
    [nH,nV] = size(M);
    a = ones(nH);

    options = odeset('RelTol',10^-8); 
    [t,x] = ode45(@ode_virus_microbe,[0 tfinal],x0,options);
    
    H = x(:,1:nH);
    V = x(:,nH+1:end);


    function xDot = ode_virus_microbe(t,x)
            H= x(1:nH);
            V = x(nH+1:end);

            Hdot = H.*(r.*(1 - a*H./K) - (M.*phi)*V);
            Vdot = V.*( (M.*phi.*beta)'*H - m);

            xDot = [Hdot;Vdot];
    end
end