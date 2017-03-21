function theta = l0_em(A, y, lambda, T)
% https://arxiv.org/pdf/1407.7508v1.pdf
% Written by Gregory Ditzler 
[n, p] = size(A);
eps_stop = .01;
eps_zero = .01;
theta = (A'*A + lambda*eye(p))\A'*y;

t = 1;
while t < T
  eta = theta;
  A_eta = repmat(eta.^2, 1, n)';
  A_eta = A_eta.*A;
  theta = (A_eta'*A + lambda*eye(p))\A_eta'*y;
  if norm(theta-eta, 2) <= eps_stop
    break
  end
end

theta(abs(theta) < eps_zero) = 0;
