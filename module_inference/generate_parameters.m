function params = generate_parameters(nH, nV, M)
%%% Generates a feasible parameter set of size nH by nV

% sampling ranges
phiMin = 10^-8; phiMax = 10^-7;
betaMin = 10; betaMax = 50;
Hmin = 10^3; Hmax = 10^4;
Vmin = 10^6; Vmax = 10^7;
K = Hmax*100;

% interaction network
if ~exist('M','var') || isempty(M)
    M = randi(2,[nH nV])-1;
end

% uniformly sample paranter ranges
phi = phiMin + (phiMax - phiMin)*rand(nH,nV);
beta = betaMin + (betaMax - betaMin)*rand(nH,nV);

% set equilibria and solve for remaining parameters
Hstar = Hmin + (Hmax-Hmin)*rand(nH,1);
Vstar = Vmin + (Vmax - Vmin)*rand(nV,1);
m = (beta.*phi.*M)'*Hstar;
r = (phi.*M)*Vstar./(1 - Hstar/K);

params = {K,r,m,M,phi,beta};

end

