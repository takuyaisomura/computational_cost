
%--------------------------------------------------------------------------------
% generate_random_coefficients.m
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

function [a1,a2,a3,a4] = generate_random_coefficients(N,sigma)

a1 = zeros(1,N^1,'single');
a2 = zeros(1,N^2,'single');
a3 = [];
a4 = [];
if (length(sigma)>=3), a3 = zeros(1,N^3,'single'); end
if (length(sigma)==4), a4 = zeros(1,N^4,'single'); end

for i1 = 1:N
 a1(1,i1) = randn(1,1,'single') * sigma(1);
 for i2 = i1+1:N
  a2(1,N*(i1-1)+i2) = randn(1,1,'single') * sigma(2);
  for i3 = i2+1:N
   if (length(sigma)>=3), a3(1,N*(N*(i1-1)+i2-1)+i3) = randn(1,1,'single') * sigma(3); end
   if (length(sigma)==4), a4(1,N*(N*(N*(i1-1)+i2-1)+i3-1)+(i3+1:N)) = randn(1,N-i3,'single') * sigma(4); end
  end
 end
end
if (length(sigma)==4), a3 = reshape(a3, [N N^2]); end
if (length(sigma)==4), a4 = reshape(a4, [N^2 N^2]); end

end

%--------------------------------------------------------------------------------
