
%--------------------------------------------------------------------------------
% cost.m
%
% This demo is included in
% Quadratic speedup of the global search using a biased crossover of two good solutions
% Takuya Isomura
%
% The MATLAB scripts are available at
% https://github.com/takuyaisomura/computational_cost
%
% Copyright (C) 2021 Takuya Isomura
% (RIKEN Center for Brain Science)
%
% 2021-10-01
%

%--------------------------------------------------------------------------------

function L = cost(x,a1,a2,a3,a4)

N    = length(x(:,1));
xx   = kron(x,ones(N,1)).*kron(ones(N,1),x);
if (isempty(a3) && isempty(a4))
 L    = a1*x + a2*xx;
elseif (isempty(a4))
 L    = a1*x + a2*xx + sum((a3*xx).*x);
else
 L    = a1*x + a2*xx + sum((a3*xx).*x) + sum((a4*xx).*xx);
end

end

%--------------------------------------------------------------------------------
