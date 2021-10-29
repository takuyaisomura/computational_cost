
%--------------------------------------------------------------------------------
% fig2b.m
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
% 2021-10-26
%

%--------------------------------------------------------------------------------
% initialisation

clear
K     = 4;         % order of cost function
N     = 10;        % state dimensionality
T     = 4000000;   % number of samples
T1    = 1000000;   % size of mini-batch
NT    = 140;       % number of sessions

seed  = 0;
rng(1000000+seed); % set seed for reproducibility

Lmean = zeros(NT,1,'single');
Lvar  = zeros(NT,1,'single');
Lmin  = zeros(NT,1,'single');
Nlist = zeros(NT,1,'single');
Tlist = zeros(NT,1,'single');
N_LM  = zeros(NT,1,'single');

%--------------------------------------------------------------------------------

for h = 1:NT
 N     = (floor((h-1)/20)+2)*5;              % state dimensionality
 T     = 4000000;                            % number of samples
 if (K <= 3 && N >= 25), T = 40000000; end   % increase T for accuracy
 if (K == 2 && N == 35), T = 400000000; end  % increase T for accuracy
 if (K == 2 && N == 40), T = 1000000000; end % increase T for accuracy
 sigma = [1 1 1 1] / sqrt(nchoosek(N,1)+nchoosek(N,2)+nchoosek(N,3)+nchoosek(N,4)); % standard derivation of coefficients
 
 % create cost function
 a1 = zeros(1,N^1,'single');
 a2 = zeros(1,N^2,'single');
 a3 = zeros(1,N^3,'single');
 a4 = zeros(1,N^4,'single');
 for i1 = 1:N
  a1(1,i1) = randn(1,1,'single') * sigma(1);
  for i2 = i1+1:N
   a2(1,N*(i1-1)+i2) = randn(1,1,'single') * sigma(2);
   for i3 = i2+1:N
    a3(1,N*(N*(i1-1)+i2-1)+i3) = randn(1,1,'single') * sigma(3);
    a4(1,N*(N*(N*(i1-1)+i2-1)+i3-1)+(i3+1:N)) = randn(1,N-i3,'single') * sigma(4);
   end
  end
 end
 a2 = reshape(a2, [N N]);
 a3 = reshape(a3, [N N^2]);
 a4 = reshape(a4, [N^2 N^2]);
 if (K == 2), a3 = []; a4 = []; end
 if (K == 3), a4 = []; end
 
%--------------------------------------------------------------------------------
 
 % global search
 L  = zeros(T,1,'single');
 for t = 0:T1:T-1
  if (rem(t,T1) == 0), fprintf('t = %d, mean = %f, var = %f, min = %f, N_LM = %d\n', t, mean(L(1:t)), var(L(1:t)), min(L(1:t)), N_LM(h,1)), end
  x      = randi([0 1],N,T1)*2 - 1;
  x1     = x;
  if (K >= 3), x2 = kron(x1,ones(N,1)).*kron(ones(N,1),x1); end
  if (K == 2), L(t+1:t+T1) = a1*x1 + sum((a2*x1).*x1);
  elseif (K == 3), L(t+1:t+T1) = a1*x1 + sum((a2*x1).*x1) + sum((a3*x2).*x1);
  elseif (K == 4), L(t+1:t+T1) = a1*x1 + sum((a2*x1).*x1) + sum((a3*x2).*x1) + sum((a4*x2).*x2); end
  L1     = L(t+1:t+T1)';
  idx    = 1:T1;
  for i = 1:N
   L1      = L1(idx);
   x1      = x1(:,idx);
   y1      = x1;
   y1(i,:) = -y1(i,:);
   if (K >= 3), y2 = kron(y1,ones(N,1)).*kron(ones(N,1),y1); end
   if (K == 2), L2 = a1*y1 + sum((a2*y1).*y1);
   elseif (K == 3), L2 = a1*y1 + sum((a2*y1).*y1) + sum((a3*y2).*y1);
   elseif (K == 4), L2 = a1*y1 + sum((a2*y1).*y1) + sum((a3*y2).*y1) + sum((a4*y2).*y2); end
   idx     = find(L1 < L2);
   if (isempty(idx) == 1), break, end
  end
  N_LM(h,1) = N_LM(h,1) + length(idx);
 end
 
%--------------------------------------------------------------------------------
 
 Nlist(h)   = N;
 Tlist(h)   = T;
 Lmean(h,1) = mean(L(1:T));
 Lvar(h,1)  = var(L(1:T));
 Lmin(h,1)  = min(L(1:T));
 
 hist(L,100)
 drawnow
 
 csvwrite(['fig2b_Lstat_K',num2str(K),'.csv'],[(1:NT)' Nlist Lmean Lvar Lmin Tlist N_LM])
end

%--------------------------------------------------------------------------------
