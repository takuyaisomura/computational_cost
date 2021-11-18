
%--------------------------------------------------------------------------------
% fig3.m
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
GRADIENT_DESCENT = 1;
K      = 4;
N      = 200;
if (GRADIENT_DESCENT == 0)
 T      = 20000;
 T1     = 1000;
 NT     = 520;
else
 T      = 100;
 NT     = 220;
end

seed   = 0;
rng(1000000+seed); % set seed for reproducibility

if (K == 2), sigma = ones(1,K) / sqrt(nchoosek(N,1)+nchoosek(N,2)); end
if (K == 3), sigma = ones(1,K) / sqrt(nchoosek(N,1)+nchoosek(N,2)+nchoosek(N,3)); end
if (K == 4), sigma = ones(1,K) / sqrt(nchoosek(N,1)+nchoosek(N,2)+nchoosek(N,3)+nchoosek(N,4)); end

Lmean  = zeros(NT,2,'single');
Lvar   = zeros(NT,2,'single');
Lmin   = zeros(NT,2,'single');
norm   = zeros(NT,1,'single');
glist  = zeros(NT,1,'single');
Lxlist = zeros(NT,T,'single');
Lylist = zeros(NT,T,'single');

%--------------------------------------------------------------------------------

for h = 1:NT
 
 % coefficients
 [a1,a2,a3,a4] = generate_random_coefficients(N,sigma);
 
 %--------------------------------------------------------------------------------
 % parents
 
 x           = randi([0 1],N,T,'single')*2 - 1;
 Lx          = zeros(T,1,'single');
 if (GRADIENT_DESCENT == 0)
  for t = 0:T1:T-1
   Lx(t+1:t+T1) = cost(x(:,t+1:t+T1),a1,a2,a3,a4);
   fprintf('.')
  end
  fprintf('\n')
 else
  [Lx,x_]     = gradient_descent(x,a1,a2,a3,a4);
 end
 
 Lmean(h,1)  = mean(Lx);
 Lvar(h,1)   = var(Lx);
 Lmin(h,1)   = min(Lx);
 Lxlist(h,:) = Lx;
 fprintf('h = %d, parents, mean = %f, var = %f, min = %f\n', h, mean(Lx), var(Lx), min(Lx))
 if (GRADIENT_DESCENT == 0)
  csvwrite(['fig3_Lstat_N', num2str(N), '.csv'],[(1:NT)' glist Lmean Lvar Lmin norm])
 else
  csvwrite(['fig3_GD_Lstat_N', num2str(N), '.csv'],[(1:NT)' glist Lmean Lvar Lmin norm])
 end
 subplot(2,1,1)
 hist(Lx,100)
 drawnow
 
 %--------------------------------------------------------------------------------
 % selection and crossover
 
 [~,idx]   = sort(Lx,'ascend');
 Py_1      = zeros(N,T,'single');
 if (GRADIENT_DESCENT == 0)
  gamma     = rem(h-1,NT/20)/50;
  for t = 1:T
   Py_1(:,t) = (x(:,idx(1))*(1-gamma)+x(:,idx(2))*gamma)/2 + 1/2;
  end
 else
  gamma     = rem(h-1,NT/20)/20;
  for t = 1:T
   Py_1(:,t) = (x_(:,idx(1))*(1-gamma)+x_(:,idx(randi([2 11])))*gamma)/2 + 1/2;
  end
 end
 
 %--------------------------------------------------------------------------------
 % offspring
 
 y           = (rand(N,T,'single') < Py_1) * 2 - 1;
 Ly          = zeros(T,1,'single');
 if (GRADIENT_DESCENT == 0)
  for t = 0:T1:T-1
   Ly(t+1:t+T1) = cost(y(:,t+1:t+T1),a1,a2,a3,a4);
   fprintf('.')
  end
  fprintf('\n')
 else
  [Ly,y_]     = gradient_descent(y,a1,a2,a3,a4);
 end
 
 Lmean(h,2)  = mean(Ly);
 Lvar(h,2)   = var(Ly);
 Lmin(h,2)   = min(Ly);
 Lylist(h,:) = Ly;
 glist(h)    = gamma;
 y1          = (2*Py_1(:,1)-1).^2;
 y2          = kron(y1,y1);
 norm(h)     = (a1.^2)*y1 + (a2.^2)*y2 + y1'*(a3.^2)*y2 + y2'*(a4.^2)*y2;
 fprintf('h = %d, offspring, mean = %f, var = %f, min = %f\n', h, mean(Ly), var(Ly), min(Ly))
 if (GRADIENT_DESCENT == 0)
  csvwrite(['fig3_Lstat_N', num2str(N), '.csv'],[(1:NT)' glist Lmean Lvar Lmin norm])
  csvwrite(['fig3_Lraw_N', num2str(N), '.csv'],[1:T 1:T; Lxlist(5:26:NT,:) Lylist(5:26:NT,:)])
 else
  csvwrite(['fig3_GD_Lstat_N', num2str(N), '.csv'],[(1:NT)' glist Lmean Lvar Lmin norm])
  csvwrite(['fig3_GD_Lraw_N', num2str(N), '.csv'],[1:T 1:T; Lxlist Lylist])
 end
 subplot(2,1,2)
 hist(Ly,100)
 drawnow
 
end

%--------------------------------------------------------------------------------
