function [Mtilde_hat, m_hat, h, W] = infer_network(t,H,V,lambda)

nT = length(t)-1;
nH = size(H,2);
nV = size(V,2);

dt = t(2:end)-t(1:end-1);
W = transpose((log(V(2:end,:))-log(V(1:end-1,:)))./dt);
h = transpose(H(2:end,:));
O = ones(1,nT);

% reconstruct infection matrix
cvx_begin
    variables Mtilde_hat(nV,nH) m_hat(nV,1)
    minimize( norm(W-Mtilde_hat*h+m_hat*O,'fro') + lambda*norm(Mtilde_hat,1) )
    subject to
    0 <= Mtilde_hat
    0 <= m_hat
cvx_end

Mtilde_hat = Mtilde_hat';

end

