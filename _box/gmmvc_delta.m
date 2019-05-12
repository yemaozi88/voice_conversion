function Y = gmmvc_delta(X,g,it)
% Y = gmmvc_delta(X,g,it)
% parameter generation considering delta-constraint
%
% NOTE
% this code is give by Miaomiao Wen, after Saito gave it to her
%
% Daisuke Saito (D3)
%

%% test
% clear all, fclose all, clc;
% % source speaker's speech
% %fin = 'J:\VoiceConversion\STRAIGHT-based\fws\scep18+delta_test\j01.feature';
% %X = loadBin(fin, 'float', 36);
% % hand motion
% fin = 'J:\H2SwithDelta\!dgvd\1\ai.dgvd';
% X = loadBin(fin, 'float', 36);
% clear fin
%                      
% % GMM object
% load('J:\H2SwithDelta\withDelta\gmm_full\jointModel_mix4_obj');
% g = obj;
% clear obj


if nargin ~= 3
    disp('usage: Y = gmmvc_delta(X,jointgmm,iter);')
    Y = 0;
    return;
end

Y = cell(1,it);

jdim = g.NDimensions;
mmix = g.NComponents;
dim = size(X,1);
T = size(X,2);

% -- delta window matrix --
% w = zeros(2*T,T); % for 1D sequence
% 
% w(1,1) = 1;
% w(2,2) = 0.5;
% for ii = 2:T-1
%     w(2*ii-1:2*ii,ii-1:ii+1) = [0 1 0;-0.5 0 0.5];
% end
% w(2*T-1,T) = 1;
% w(2*T,T-1) = -0.5;

w = sparse(dim*T,dim*T/2);
for ii = 1:T
    w((ii-1)*dim+1:(2*ii-1)*dim/2,(ii-1)*dim/2+1:ii*dim/2) = eye(dim/2);
end
for ii = 1:T-1
    w((2*ii-1)*dim/2+1:ii*dim,ii*dim/2+1:(ii+1)*dim/2) = 0.5 * eye(dim/2);
end
for ii = 2:T
    w((2*ii-1)*dim/2+1:ii*dim,(ii-2)*dim/2+1:(ii-1)*dim/2) = -0.5 * eye(dim/2);
end


xmu = g.mu(:,1:dim); % mmix by d
xsigma = zeros(1,dim,mmix);
weightsigma = zeros(1,dim,mmix);
ymu = g.mu(:,dim+1:jdim);
ysigma = zeros(1,dim,mmix);
jsigma = zeros(1,dim,mmix);

if dim*2 ~= jdim
    error('Dim of X should be half of jointmodel: X(%d) g1(%d)',dim,jdim);
end

jgsigma = repmat(eye(dim,dim),[2 2 mmix]) .* g.Sigma; % masking

for ii = 1:mmix
    xsigma(:,:,ii) = diag(g.Sigma(1:dim,1:dim,ii))';
    ysigma(:,:,ii) = diag(g.Sigma(dim+1:jdim,dim+1:jdim,ii))';
    weightsigma(:,:,ii) = diag(g.Sigma(dim+1:jdim,1:dim,ii))' ./ xsigma(:,:,ii);
    jsigma(:,:,ii) = ysigma(:,:,ii) - diag(g.Sigma(dim+1:jdim,1:dim,ii))' ...
        ./ xsigma(:,:,ii) .* diag(g.Sigma(1:dim,dim+1:jdim,ii))';
end

xg = gmdistribution(xmu,xsigma,g.PComponents);
yg = gmdistribution(ymu,ysigma,g.PComponents);
jg = gmdistribution(g.mu,jgsigma,g.PComponents);

gamm = posterior(xg,X'); % T by mmix

tmpdiff = repmat(X,[1,1,mmix]) - repmat(permute(xmu,[2 3 1]),[1,T,1]); % d by T by mmix
outmean1 = repmat(permute(weightsigma,[2,1,3]),[1,T,1]) .* tmpdiff ...
        + repmat(permute(ymu,[2 3 1]),[1,T,1]); % d by T by mmix
outv1 = repmat(permute(1 ./ jsigma,[3 1 2]),[1,T,1]);  %mmix by T by d

outmean2 =  permute(outmean1,[3 2 1]) .* outv1; % mmix by T by d 


fprintf('it\tlog-likelihood\n\n');
for ii = 1:it
    Y{ii} = zeros(dim/2,T);
    % mean
    tmp1 = repmat(gamm',[1,1,dim]) .* outmean2 ; % mmix by T by d
    y1 = permute(sum(tmp1,1),[3 2 1]); % 1 by T by d -> d by T
    % cov
    tmpv1 = repmat(gamm',[1,1,dim]) .* outv1 ; %mmix by T by d
    v1 = permute(sum(tmpv1,1),[3 2 1]); % 1 by T by d -> d by T
    
%     for jj = 1:dim/2
%         ytmp = (w'*diag(reshape(v1([jj,jj+dim/2],:),2*T,1))*w)\(w'*reshape(y1([jj,jj+dim/2],:),2*T,1));
%         Y{ii}(jj,:) = ytmp';
%     end
    ytmp = (w'*sparse(1:dim*T,1:dim*T,reshape(v1,dim*T,1))*w)\(w'*reshape(y1,dim*T,1));
    
    Y{ii} = reshape(ytmp,dim/2,T);    
    
    gamm = posterior(jg,[X;adddelta(Y{ii})]'); % T by mmix
    fprintf('%d\t%e\t\n',ii,sum(log(pdf(jg,[X;adddelta(Y{ii})]'))) - sum(log(pdf(xg,X'))));
end
