function Y = gmmvc_(X,g)
% Y = gmmvc_(X, g)
% calculates converted time series using conversion model
% g is trained with the same parallel data to one for gmmvc.m
% however augumented vector z = [Y, X], instead of [X, Y]
%
% X: input (without power)
% g: jointgmm (conversion model)
% Y: output (without power)
%
% HISTORY
% 2011/02/02 modified spkmodel_vc2.m for combinedModel.m
%
% AUTHOR
% Aki Kunikoshi (D2)
% yemaozi88@gmail.com
%

if nargin ~= 2
    disp('usage: Y = gmmvc(X,jointgmm);')
    Y = 0;
    return;
end

Y = zeros(size(X));

jdim = g.NDimensions;
k = g.NComponents;
dim = size(X,1);
N = size(X,2);

if dim*2 ~= jdim
    error('Dim or X should be half of jointmodel: X(%d) g(%d)',dim,jdim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% switched S2Hmodel into H2Smodel %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gmu = zeros(k, jdim);
gmu(:, 1:dim) = g.mu(:, dim+1:jdim);
gmu(:, dim+1:jdim) = g.mu(:, 1:dim);

gSigma = zeros(jdim, jdim, k);
for ii = 1:k
    gSigma(1:dim, 1:dim, ii)           = g.Sigma(dim+1:jdim, dim+1:jdim, ii);
    gSigma(dim+1:jdim, dim+1:jdim, ii) = g.Sigma(1:dim, 1:dim, ii);
    gSigma(1:dim, dim+1:jdim, ii)      = g.Sigma(dim+1:jdim, 1:dim, ii)';
    gSigma(dim+1:jdim, 1:dim, ii)      = g.Sigma(1:dim, dim+1:jdim, ii)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% modified %%%%%%%%%%%%%%%%%%
%xmu = g.mu(:,1:dim);
%ymu = g.mu(:,dim+1:jdim);
xmu = gmu(:,1:dim);
ymu = gmu(:,dim+1:jdim);

xsigma = zeros(1,dim,k);
weightsigma = zeros(1,dim,k);
xp = g.PComponents;

for ii = 1:k
    %%%%%%%%%%%%%%%%% modified %%%%%%%%%%%%%%%%%%
    %xsigma(:,:,ii) = diag(g.Sigma(1:dim,1:dim,ii))';
    %weightsigma(:,:,ii) = diag(g.Sigma(dim+1:jdim,1:dim,ii))' ./ xsigma(:,:,ii);
    xsigma(:,:,ii) = diag(gSigma(1:dim,1:dim,ii))';
    weightsigma(:,:,ii) = diag(gSigma(dim+1:jdim,1:dim,ii))' ./ xsigma(:,:,ii);

%    ysigma(:,:,ii) = diag(g.Sigma(dim+1:jdim,dim+1:jdim,ii))';
%    xysigma(:,:,ii) = diag(g.Sigma(1:dim,dim+1:jdim,ii))'; 
end

xg = gmdistribution(xmu,xsigma,xp);
px = posterior(xg,X'); % N by k


tmpdiff = repmat(X,[1,1,k]) - repmat(permute(xmu,[2 3 1]),[1,N,1]); % d by N by k
outmean = repmat(permute(weightsigma,[2,1,3]),[1,N,1]) .* tmpdiff ...
    + repmat(permute(ymu,[2 3 1]),[1,N,1]); % d by N by k

tmp = repmat(px',[1,1,dim]) .* permute(outmean,[3 2 1]); % k by N by d
Y = permute(sum(tmp,1),[3 2 1]); % 1 by N by d -> d by N