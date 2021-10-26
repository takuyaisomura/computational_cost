
%--------------------------------------------------------------------------------
% fig1.m
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
% initialisation

clear
K        = 4;      % degree of products
NT       = 240;    % number of trials
NTex     = 80;     % number of trials for brute force search

seed     = 0;
rng(1000000+seed); % set seed for reproducibility

Lmean    = zeros(NT,1,'single');
Lvar     = zeros(NT,1,'single');
Lmin     = zeros(NT,1,'single');
nlist    = zeros(NT,1,'single');
Ngdmean  = zeros(NT,1,'single');
Lmean_ex = zeros(NTex,1,'single');
Lvar_ex  = zeros(NTex,1,'single');
Lmin_ex  = zeros(NTex,1,'single');
nlist_ex = zeros(NTex,1,'single');

%--------------------------------------------------------------------------------

for h = 1:NT
 
 N      = floor((h-1)/20+1)*5;
 T      = 2^N;
 Nterms = 0;
 for i = 1:K, Nterms = Nterms + nchoosek(N,i); end
 sigma  = [ones(1,K)/sqrt(Nterms) zeros(1,4-K)];
 
%--------------------------------------------------------------------------------
 
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
 a3 = reshape(a3, [N N^2]);
 a4 = reshape(a4, [N^2 N^2]);
 if (K == 2), a3 = []; a4 = []; end
 if (K == 3), a4 = []; end
 
%--------------------------------------------------------------------------------
 
 if (h <= 80)
  if (rem(h,20) == 1)
   x  = cast([-1 1],'single');
   for i = 2:N, x = [x x; -ones(1,2^(i-1)) ones(1,2^(i-1))]; end
  end
  L             = cost(x,a1,a2,a3,a4);
  nlist_ex(h,1) = N;
  Lmean_ex(h,1) = mean(L);
  Lvar_ex(h,1)  = var(L);
  Lmin_ex(h,1)  = min(L);
  if (h == 80)
   [~,idx] = min(L);
   csvwrite(['Lstat_fig1_x_L', num2str(K),'.csv'],[x(:,1:100)' L(1:100)'; x(:,(idx-99):(idx+100))' L((idx-99):(idx+100))'; x(:,(T-99):T)' L((T-99):T)'])
   csvwrite(['Lstat_fig1_Lhist', num2str(K),'.csv'],[(-10:0.1:(10-0.1)); histcounts(L,-10:0.1:10)])
   csvwrite(['Lstat_fig1_K', num2str(K),'ex.csv'],[(1:NTex)' nlist_ex Lmean_ex Lvar_ex Lmin_ex])
  end
 end
 
 T            = 10000;
 x0           = cast(randi([0 1],N,T)*2-1,'single');
 [L,x_,Ngd]   = gradient_descent(x0,a1,a2,a3,a4);
 nlist(h,1)   = N;
 Lmean(h,1)   = mean(L);
 Lvar(h,1)    = var(L);
 Lmin(h,1)    = min(L);
 Ngdmean(h,1) = mean(Ngd);
 fprintf('h = %d, N = %d, mean = %f, var = %f, min = %f, Ngd = %f\n', h, N, Lmean(h,1), Lvar(h,1), Lmin(h,1), Ngdmean(h,1))
 
 hist(L,100)
 title(['h=',num2str(h),', N=',num2str(N)])
 drawnow
 
 csvwrite(['Lstat_fig1_K', num2str(K),'.csv'],[(1:NT)' nlist Lmean Lvar Lmin Ngdmean])
 
end

%--------------------------------------------------------------------------------

