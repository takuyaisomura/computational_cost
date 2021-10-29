
%--------------------------------------------------------------------------------
% fig2a.m
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
% 2021-10-26
%

%--------------------------------------------------------------------------------
% initialisation

clear
N    = 1002;
num  = 100;
seed = 0;
rng(1000000+seed); % set seed for reproducibility

%--------------------------------------------------------------------------------

X    = zeros(N,N);
pos  = rand(num,2)*N;
pos  = min(max(ceil(pos),2),N-1);
for k = 1:num
 for i = 1:N
  for j = 1:N, X(i,j) = X(i,j) + exp(-((i-pos(k,1))^2+(j-pos(k,2))^2)/2/50^2); end
 end
end
X = (X - min(min(X)))/(max(max(X)) - min(min(X)));
X = 1 - X;

%--------------------------------------------------------------------------------

img = zeros(N,N,3);
img(:,:,1) = max(min(X*3-2,1),0);
img(:,:,2) = max(min(X*3-1,1),0);
img(:,:,3) = 1;
image(img)
drawnow
imwrite(img(2:N-1,2:N-1,:), ['fig2a_cost_function.png'])
dlmwrite(['fig2a_cost_function.csv'],X(2:N-1,2:N-1),'precision',6)

%--------------------------------------------------------------------------------

plm = [];
for i = 1:N
 for j = 1:N
  if (i<2||j<2||i>N-1||j>N-1), continue, end
  Xmin = min([X(i,j),X(i-1,j),X(i+1,j),X(i,j-1),X(i,j+1)]);
  if (Xmin == X(i,j)), plm = [plm; i, j]; end
 end
end
NLM  = length(plm(:,1));
perm = randperm(NLM);
plm  = plm(perm,:);
Llm  = zeros(NLM,1);
for k = 1:NLM, Llm(k) = X(plm(k,1),plm(k,2)); end
dlmwrite(['fig2a_local_minima.csv'],[(1:NLM)',plm-1,Llm],'precision',6)

%--------------------------------------------------------------------------------

Y    = zeros(N,N);
for k = 1:NLM, Y(plm(k,1),plm(k,2)) = k; end
for D = 1:N/2
 for k = 1:NLM
  for l = 1:4
   for d = [0:-1:-D,1:D]
    if (l == 1), i = plm(k,1)+d; j = plm(k,2)+D;
    elseif (l == 2), i = plm(k,1)+d; j = plm(k,2)-D;
    elseif (l == 3), i = plm(k,1)+D; j = plm(k,2)+d;
    elseif (l == 4), i = plm(k,1)-D; j = plm(k,2)+d; end
    if (i<2||j<2||i>N-1||j>N-1||Y(i,j)>0), continue, end
    Xmin = min([X(i,j),X(i-1,j),X(i+1,j),X(i,j-1),X(i,j+1)]);
    if (X(i-1,j) == Xmin && Y(i-1,j) > 0), Y(i,j) = Y(i-1,j);
    elseif (X(i+1,j) == Xmin && Y(i+1,j) > 0), Y(i,j) = Y(i+1,j);
    elseif (X(i,j-1) == Xmin && Y(i,j-1) > 0), Y(i,j) = Y(i,j-1);
    elseif (X(i,j+1) == Xmin && Y(i,j+1) > 0), Y(i,j) = Y(i,j+1); end
   end
  end
 end
 img2 = hsv2rgb(Y/NLM,ones(N,N),ones(N,N));
 img2(:,:,1) = img2(:,:,1) .* (Y>0);
 img2(:,:,2) = img2(:,:,2) .* (Y>0);
 img2(:,:,3) = img2(:,:,3) .* (Y>0);
 image(img2)
 drawnow
end

%--------------------------------------------------------------------------------

for D = 1:N
 for l = 1:4
  for d = 1:N
   if (l == 1), i = d; j = D;
   elseif (l == 2), i = d; j = N-D+1;
   elseif (l == 3), i = D; j = d;
   elseif (l == 4), i = N-D+1; j = d; end
   if (i<2||j<2||i>N-1||j>N-1||Y(i,j)>0), continue, end
   Xmin = min([X(i,j),X(i-1,j),X(i+1,j),X(i,j-1),X(i,j+1)]);
   if (X(i-1,j) == Xmin && Y(i-1,j) > 0), Y(i,j) = Y(i-1,j);
   elseif (X(i+1,j) == Xmin && Y(i+1,j) > 0), Y(i,j) = Y(i+1,j);
   elseif (X(i,j-1) == Xmin && Y(i,j-1) > 0), Y(i,j) = Y(i,j-1);
   elseif (X(i,j+1) == Xmin && Y(i,j+1) > 0), Y(i,j) = Y(i,j+1); end
  end
 end
 img2 = hsv2rgb(Y/NLM,ones(N,N),ones(N,N));
 img2(:,:,1) = img2(:,:,1) .* (Y>0);
 img2(:,:,2) = img2(:,:,2) .* (Y>0);
 img2(:,:,3) = img2(:,:,3) .* (Y>0);
 image(img2)
 drawnow
end
imwrite(img2(2:N-1,2:N-1,:), ['fig2a_basins.png'])

%--------------------------------------------------------------------------------

bound = [];
for i = 3:N-2
 for j = 3:N-2
  if (Y(i,j)~=Y(i-1,j)||Y(i,j)~=Y(i+1,j)||Y(i,j)~=Y(i,j-1)||Y(i,j)~=Y(i,j+1)), bound = [bound; i j]; end
 end
end
img3  = insertShape(img2,'FilledCircle',[plm(:,2) plm(:,1) ones(NLM,1)*8],'Color','black','Opacity',1);
img3  = insertShape(img3,'FilledCircle',[bound(:,2) bound(:,1) ones(length(bound(:,1)),1)*2],'Color','black','Opacity',1);
image(img3);
drawnow
imwrite(img3(2:N-1,2:N-1,:), ['fig2a_basins2.png'])

%--------------------------------------------------------------------------------

Z = zeros(N,N);
for i = 2:N-1
 for j = 2:N-1
  Z(i,j) = X(plm(Y(i,j),1),plm(Y(i,j),2));
 end
end
img4 = zeros(N,N,3);
img4(:,:,1) = max(min(Z*3-2,1),0);
img4(:,:,2) = max(min(Z*3-1,1),0);
img4(:,:,3) = 1;
imwrite(img4(2:N-1,2:N-1,:), ['fig2a_cost_function2.png'])
dlmwrite(['fig2a_cost_function2.csv'],Z(2:N-1,2:N-1),'precision',6)

%--------------------------------------------------------------------------------

img5  = insertShape(img4,'FilledCircle',[bound(:,2) bound(:,1) ones(length(bound(:,1)),1)*1],'Color','black','Opacity',0.2);
imwrite(img5(2:N-1,2:N-1,:), ['fig2a_cost_function3.png'])

fig = figure();
subplot(3,2,3)
surf(X(2:N-1,2:N-1),img(2:N-1,2:N-1,:),'EdgeColor','none')
camproj('perspective')
ax            = gca;
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.ZTickLabel = [];
ax.LineWidth  = 2;
subplot(3,2,4)
surf(Z(2:N-1,2:N-1),img5(2:N-1,2:N-1,:),'EdgeColor','none')
camproj('perspective')
ax            = gca;
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.ZTickLabel = [];
ax.LineWidth  = 2;
set(fig, 'PaperPosition', [0,0,20,28])
print(fig, 'fig2c_cost_function.png', '-dpng', '-r600');

%--------------------------------------------------------------------------------

