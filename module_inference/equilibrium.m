function x = equilibrium(params)
% equilibrium(M,r,a,K,phi,beta,m)

[K,r,m,M,phi,beta] = params{:};
[nH,nV] = size(M);
a = ones(nH);
cMatrix = [ a.*repmat(r,1,nH)/K, M.*phi; (M.*phi.*beta)', zeros(nV)];
x = cMatrix\[r;m];

end