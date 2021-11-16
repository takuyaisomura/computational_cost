
%--------------------------------------------------------------------------------
% fig4.m
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
% initialisation

clear
tic
N      = 500;                  % state dimensionality
T      = 500;                  % number of iterations
NT     = 50;                   % number of sessions

seed   = 0;
rng(1000000+seed);             % set seed for reproducibility

Lxlist = zeros(NT,T,'single'); % list of final costs in 1st generation
Lylist = zeros(NT,T,'single'); % list of final costs in 2nd generation

%--------------------------------------------------------------------------------

% run simulation
for h = 1:NT
 pos = rand(2,N,'single');      % positions of nodes
 A   = zeros(N,N,'single');     % cost = distance
 for i = 1:N, A(i,:) = sqrt((pos(1,i)-pos(1,:)).^2 + (pos(2,i)-pos(2,:)).^2); end
 A0  = A;
 A   = (A - mean(A(A~=0))) / std(A(A~=0)) / sqrt(N);
 A   = A - diag(diag(A));
 
 %--------------------------------------------------------------------------------
 
 % parents
 fprintf(1,'[h = %d, parents]\n', h)
 x0  = zeros(N,T,'single');     % initial states
 x   = zeros(N,T,'single');     % final states
 Lx0 = zeros(1,T,'single');     % initial costs
 Lx  = zeros(1,T,'single');     % final costs
 
 for t = 1:T
  fprintf(1,'t = %d (%.0f sec) ', t, toc)
  x0(:,t) = [1 randperm(N-1)+1]';                % initial state
  Lx0(t)  = sum(A(x0(:,t)+N*(x0([2:N 1],t)-1))); % initial cost
  xt      = x0(:,t);                             % state
  Diff    = tsp_gradient(xt,A,[],1,N);           % gradient
  for rep = 1:N^2
   if (rem(rep,100) == 0), fprintf(1,'.'), end
   [d_min,idx_min] = min(Diff(1:N^2));           % detect steepest descent direction
   if (d_min >= 0), break, end                   % if the minimum gradient is non-negative, break the loop
   i_min = rem(idx_min-1,N)+1;                   % start point
   j_min = floor((idx_min-1)/N)+1;               % end point
   xt(i_min:j_min-1) = flip(xt(i_min:j_min-1));  % change the route
   Diff = tsp_gradient(xt,A,Diff,i_min,j_min);   % update the gradient
  end
  fprintf(1,'\n')
  x(:,t)  = xt;                                  % final state
  Lx(t)   = sum(A(x(:,t)+N*(x([2:N 1],t)-1)));   % final cost
 end
 
 subplot(2,2,1)
 hist(Lx0,100)
 title('1st generation initial costs')
 subplot(2,2,2)
 hist(Lx,100)
 title('1st generation final costs')
 drawnow
 
 %--------------------------------------------------------------------------------
 
 % selection and crossover
 fprintf(1,'[h = %d, selection and crossover]\n', h)
 y0  = zeros(N,T,'single');     % initial states
 y   = zeros(N,T,'single');     % final states
 Ly0 = zeros(1,T,'single');     % initial costs
 Ly  = zeros(1,T,'single');     % final costs
 
 Ntop    = 10;
 [~,idx] = sort(Lx,'ascend');
 eps     = 0.05;                % crossover rate
 
 temp          = zeros(N,T);
 num_mismatch0 = zeros(T,1);
 num_mismatch1 = zeros(T,1);
 for t = 1:T
  i1    = randi([1 Ntop]);
  i2    = randi([1 Ntop]);
  while i1 == i2, i2 = randi([1 Ntop]); end
  Map1  = zeros(N,N);
  Map2  = zeros(N,N);
  Map3  = zeros(N,N);
  Map1(x(:,idx(i1))+N*(x([2:N 1],idx(i1))-1)) = 1;
  Map2(x(:,idx(i2))+N*(x([2:N 1],idx(i2))-1)) = 1;
  rnd_e = (rand(N,1) < eps) * 1;
  for i = 1:N
   if (rnd_e(i) == 0), Map3(i,:) = Map1(i,:);
   else, Map3(i,:) = Map2(i,:); end
  end
  rnd_i = randperm(N)';
  for i = 1:N
   j = (1:N)*Map3(rnd_i(i),:)';
   if (sum(Map3(:,j)) > 1)
    if (rnd_e(rnd_i(i)) == 0), Map3(rnd_i(i),:) = Map2(rnd_i(i),:);
    else, Map3(rnd_i(i),:) = Map1(rnd_i(i),:); end
   end
  end
  rnd_j = randperm(N)';
  y0(1,t) = 1;
  temp(:,t) = 0;
  temp(1,t) = 1;
  for i = 2:N
   y0(i,t) = (1:N)*Map3(y0(i-1,t),:)';
   if (temp(y0(i,t),t) > 0)
    for k = 1:N
     if (sum(Map3(:,rnd_j(k))) == 0 && temp(rnd_j(k),t) == 0), break, end
    end
    if (temp(rnd_j(k),t) > 0)
     for k = 1:N
      if (temp(rnd_j(k),t) == 0), break, end
     end
    end
    y0(i,t) = rnd_j(k);
    Map3(y0(i-1,t),:) = 0;
    Map3(y0(i-1,t),y0(i,t)) = 1;
   end
   temp(y0(i,t),t) = temp(y0(i,t),t) + 1;
  end
  Map3(y0(N,t),:) = 0;
  Map3(y0(N,t),y0(1,t)) = 1;
  num_mismatch0(t) = sum(sum(Map1 ~= Map2));
  num_mismatch1(t) = sum(sum(Map1 ~= Map3));
  fprintf(1,'.')
 end
 fprintf(1,'\n')
 
 %--------------------------------------------------------------------------------
 
 % offspring
 fprintf(1,'[h = %d, offspring]\n', h)
 for t = 1:T
  fprintf(1,'t = %d (%.0f sec) ', t, toc)
  Ly0(t)  = sum(A(y0(:,t)+N*(y0([2:N 1],t)-1))); % initial cost
  yt      = y0(:,t);                             % state
  Diff    = tsp_gradient(yt,A,[],1,N);           % gradient
  for rep = 1:N^2
   if (rem(rep,100) == 0), fprintf(1,'.'), end
   [d_min,idx_min] = min(Diff(1:N^2));           % detect steepest descent direction
   if (d_min >= 0), break, end                   % if the minimum gradient is non-negative, break the loop
   i_min = rem(idx_min-1,N)+1;                   % start point
   j_min = floor((idx_min-1)/N)+1;               % end point
   yt(i_min:j_min-1) = flip(yt(i_min:j_min-1));  % change the route
   Diff = tsp_gradient(yt,A,Diff,i_min,j_min);   % update the gradient
  end
  fprintf(1,'\n')
  y(:,t)  = yt;                                  % final state
  Ly(t)   = sum(A(y(:,t)+N*(y([2:N 1],t)-1)));   % final cost
 end
 
 subplot(2,2,3)
 hist(Ly0,100)
 title('2nd generation initial costs')
 subplot(2,2,4)
 hist(Ly,100)
 title('2nd generation final costs')
 drawnow
 
 %--------------------------------------------------------------------------------
 
 Lxlist(h,:) = Lx;
 Lylist(h,:) = Ly;
 csvwrite(['TSP_Lraw_N', num2str(N), '.csv'],[1:T 1:T; Lxlist Lylist])
 
 if (h == 1)
  [~,idx]  = sort(Lx,'ascend');
  [~,idx2] = sort(Ly,'ascend');
  csvwrite(['TSP_Route_N', num2str(N), '.csv'],[(1:N)' pos(:,x0(:,1))' pos(:,x(:,idx(1)))' pos(:,x(:,idx(2)))' pos(:,y(:,idx2(1)))'])
 end
 fprintf(1,'----------------------------------------\n')
end

%--------------------------------------------------------------------------------

