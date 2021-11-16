
%--------------------------------------------------------------------------------
% tsp_gradient.m
%
% This demo is included in
% Quadratic speedup of global search using a biased crossover of two good solutions
% Takuya Isomura
%
% The MATLAB scripts are available at
% https://github.com/takuyaisomura/computational_cost
%
% Copyright (C) 2021 Takuya Isomura
% (RIKEN Center for Brain Science)
%
% 2021-08-12
%

%--------------------------------------------------------------------------------

function Diff = tsp_gradient(xt,A,Diff0,i_,j_)
N      = length(xt);
if (isempty(Diff0) == 1), Diff = zeros(N,N);
else, Diff = Diff0; end
xt_1   = xt([N 1:N-1]);
for ii = i_:j_
 i    = [ii*ones(1,length(ii+1:N)),1:i_-1];
 j    = [ii+1:N,ii*ones(1,length(1:i_-1))];
 idx1 = [xt_1(i)+N*(xt(i)-1), xt_1(j)+N*(xt(j)-1)]';
 idx2 = [xt_1(i)+N*(xt_1(j)-1), xt(i)+N*(xt(j)-1)]';
 Diff(i+N*(j-1)) = sum(A(idx2)) - sum(A(idx1));
end
end

