
%--------------------------------------------------------------------------------
% gradient_descent.m
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

function [L, xt, Ngd] = gradient_descent(x,a1,a2,a3,a4)

N    = length(x(:,1));
T    = length(x(1,:));
Ngd  = zeros(T,1,'single');
xt   = x;

%--------------------------------------------------------------------------------

% gradient descent algorithm
for rep = 1:N
 
 % compute current cost
 x1   = xt;
 x2   = kron(x1,ones(N,1)).*kron(ones(N,1),x1);
 a1x1 = a1 * x1;
 a2x2 = a2 * x2;
 if (isempty(a3) && isempty(a4))
  L    = a1x1 + a2x2;
 elseif (isempty(a4))
  a3x2 = a3 * x2;
  L    = a1x1 + a2x2 + sum(x1.*a3x2);
 else
  a3x2 = a3 * x2;
  a4x2 = a4 * x2;
  L    = a1x1 + a2x2 + sum(x1.*a3x2) + sum(x2.*a4x2);
 end
 
 % compute gradient
 DL   = zeros(T,N,'single');
 for i = 1:N
  Dx1      = x1;
  Dx1(i,:) = -Dx1(i,:);
  Dx2      = kron(Dx1,ones(N,1)).*kron(ones(N,1),Dx1);
  Dx2_x2   = Dx2 - x2;
  idx      = find(Dx2_x2(:,1) ~= 0);
  a1Dx1    = a1 * Dx1;
  a2Dx2    = a2 * Dx2;
  if (isempty(a3) && isempty(a4))
   DL(:,i)  = a1Dx1 + a2Dx2;
  elseif (isempty(a4))
   a3Dx2    = a3(:,idx) * Dx2_x2(idx,:) + a3x2;
   DL(:,i)  = a1Dx1 + a2Dx2 + sum(Dx1.*a3Dx2);
  else
   a3Dx2    = a3(:,idx) * Dx2_x2(idx,:) + a3x2;
   a4Dx2    = a4(:,idx) * Dx2_x2(idx,:) + a4x2;
   DL(:,i)  = a1Dx1 + a2Dx2 + sum(Dx1.*a3Dx2) + sum(Dx2.*a4Dx2);
  end
 end
 
 % comparison and state update
 flag_continue = 0;
 for t = 1:T
  [Lmin,idx] = min(DL(t,:));
  if (Lmin < L(t))
   xt(idx,t) = -xt(idx,t);
   flag_continue = 1;
   Ngd(t) = rep;
  end
 end
 fprintf(1,'.')
 if (flag_continue == 0), break, end
 
end
fprintf(1,'\n')

end

%--------------------------------------------------------------------------------
