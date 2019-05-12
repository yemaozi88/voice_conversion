function Y = gmmMappingWithDelta(X, obj)
% Y = gmmMappingWithDelta(X, obj)
% estimates Y using obj trained by parallel data of [X, Y] which include
% delta features (static part and delta part are separated)
%
% INPUT
% X: input (without power)
% g: jointgmm (conversion model, full covariance)
% Y: output (without power)
%
% NOTE
% - refered to...
% Tomoki Toda, Alan W Black, Keiichi Tokuda:
% "Voice Conversion Based on Maximum Likelihood Estimation of Spectral Parameter Trajectory"
%
% HISTORY
% 2011/08/18 the bug that overwriting sigmaXX with sigmaYY in loading
% covariance (diag) part is fixed
% 2011/08/11 functionized
%
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%


%% test
% clear all, fclose all, clc;
% % source speaker's speech
% fin = 'J:\VoiceConversion\STRAIGHT-based\fws\scep18+delta_test\j01.feature';
% X = loadBin(fin, 'float', 36);
% clear fin
% 
% % GMM object
% load('C:\research\VoiceConversion\STRAIGHT-based\fws2mht_withDelta_reduced1of5\gmm_full\jointModel_mix8_obj');


%% load GMM data
jDim = obj.NDimensions; % the number of dimensions of joint vector
gNum = obj.NComponents; % the number of mixtures: k
sDim = size(X, 1); % the number of dimensions of source data: dim
sNum = size(X, 2); % the number of frames in source data: N

dDim = sDim/2; % the number of dimensions of static/delta part 
if floor(dDim)*2 ~= sDim
    error('spectral part and delta part should have the same length');
end
x_idx = [1:sDim];      % source
y_idx = [sDim+1:jDim]; % target
s_idx = [1:dDim];      % static part
d_idx = [dDim+1:sDim]; % delta part

% mean
muX = obj.mu(:, x_idx); % gNum x sDim
muY = obj.mu(:, y_idx); % gNum x tDim (jDim - sDim)
p = obj.PComponents; % NComponents x 1, proportions

% covariance (full)
% sigmaXX = zeros(sDim, sDim, gNum);
% sigmaXY = zeros(sDim, sDim, gNum);
% sigmaYX = zeros(sDim, sDim, gNum);
% sigmaYY = zeros(sDim, sDim, gNum);
% for ii = 1:gNum
%     sigmaXX(:, :, ii) = obj.Sigma(x_idx, x_idx, ii);
%     sigmaXY(:, :, ii) = obj.Sigma(x_idx, y_idx, ii);
%     sigmaYX(:, :, ii) = obj.Sigma(y_idx, x_idx, ii);
%     sigmaYY(:, :, ii) = obj.Sigma(y_idx, y_idx, ii);
% end

% covariance (diag)
sigmaXX = zeros(1, sDim, gNum);
sigmaXY = zeros(1, sDim, gNum);
sigmaYX = zeros(1, sDim, gNum);
sigmaYY = zeros(1, sDim, gNum);
for ii = 1:gNum
    sigmaXX(:, :, ii) = diag(obj.Sigma(x_idx, x_idx, ii))';
    sigmaXY(:, :, ii) = diag(obj.Sigma(x_idx, y_idx, ii))';
    sigmaYX(:, :, ii) = diag(obj.Sigma(y_idx, x_idx, ii))';
    sigmaYY(:, :, ii) = diag(obj.Sigma(y_idx, y_idx, ii))';
%    sigmaXX(:, :, ii) = diag(obj.Sigma(1, x_idx, ii))';
%    sigmaXY(:, :, ii) = diag(obj.Sigma(1, y_idx, ii))';
%    sigmaYX(:, :, ii) = diag(obj.Sigma(1, x_idx, ii))';
%    sigmaYY(:, :, ii) = diag(obj.Sigma(1, y_idx, ii))';
end
clear x_idx y_idx obj


%% E and D matrice
% Y = zeros(size(X));
objX = gmdistribution(muX, sigmaXX, p);

for ii = 1:sNum
    x = X(:, ii);

    % choose the most likely distribution
    pX = posterior(objX, x'); % X should be fNum x dim
    [maxValue, jj] = max(pX);
    clear pX maxValue
    
    sigmaXX_ = squeeze(sigmaXX(:, :, jj))';
    sigmaXY_ = squeeze(sigmaXY(:, :, jj))';
    sigmaYX_ = squeeze(sigmaYX(:, :, jj))';
    sigmaYY_ = squeeze(sigmaYY(:, :, jj))';
    
    invSigmaXX_ = 1 ./ sigmaXX_;
    %E = muY(jj, :)' + sigmaYX_ * inv(sigmaXX_) * (x - muX(jj, :)');
    %D = sigmaYY_ - sigmaYX_ * inv(sigmaXX_) * sigmaXY_;
    E = muY(jj, :)' + sigmaYX_ .* invSigmaXX_ .* (x - muX(jj, :)');
    D = sigmaYY_ - sigmaYX_ .* invSigmaXX_ .* sigmaXY_;

    E_s(:, ii) = E(s_idx); % mean of spectral part
    E_d(:, ii) = E(d_idx); % mean of delta spectral part
    D_s(:, ii) = D(s_idx); % sigma of spectral part
    D_d(:, ii) = D(d_idx); % sigma of delta spectral part
    
    clear sigmaXX_ sigmaXY_ sigmaYX_ sigmaYY_
    clear E D
end %sNum
clear x ii jj s_idx d_idx
clear muX muY p objX
clear sigmaXX sigmaXY sigmaYX sigmaYY

Evec_s = reshape(E_s, dDim * sNum, 1);
Evec_d = reshape(E_d, dDim * sNum, 1);
Dvec_s = reshape(D_s, dDim * sNum, 1);
Dvec_d = reshape(D_d, dDim * sNum, 1);
clear E_s E_d D_s D_d;


%% W matrix
row = sNum * dDim;
col = (sNum + 2) * dDim;

% weight: static: 1, delta: 0.5

% static
W_s = spalloc(row, col, col); % sparce marix which has 'col' non-zero factors maximum
for h = 1: sNum*dDim
    W_s(h, h + dDim) = 1;
end
clear h

% delta
W_d = spalloc(row, col, col*2);    
for t = 1:sNum
  for p = 1:dDim      
     W_d((t-1)*dDim+p, (t-1)*dDim+p) = -0.5;
     W_d((t-1)*dDim+p, (t+1)*dDim+p) = 0.5;              
  end     
end 
clear t p
clear row col


%% mapping
Dm_s = spalloc(dDim*sNum, dDim*sNum, dDim*sNum);
Dm_d = spalloc(dDim*sNum, dDim*sNum, dDim*sNum);

for mDim = 1 : dDim*sNum
     Dm_s(mDim, mDim) = Dvec_s(mDim);
     Dm_d(mDim, mDim) = Dvec_d(mDim);
end
clear mDim

invDm_s = inv(Dm_s);
invDm_d = inv(Dm_d);
clear Dvec_s Dvec_d;

% eq. 39
traW_invDm_s = W_s' * invDm_s;
traW_invDm_d = W_d' * invDm_d;

A_s = traW_invDm_s * W_s;
A_d = traW_invDm_d * W_d;

b_s = traW_invDm_s * Evec_s;
b_d = traW_invDm_d * Evec_d;

Y = (A_s + A_d) \ (b_s + b_d); 

clear Evec_s Evec_d W_s W_d
clear Dm_s Dm_d invDm_s invDm_d traW_invDm_s traW_invDm_d
clear A_s A_d b_s b_d 

% reshape
Y = reshape(Y, dDim, length(Y)/dDim);
clear dDim sDim jDim gNum sNum 


%% check
% % scep
% fsrcScep = 'J:\VoiceConversion\STRAIGHT-based\fws\scep18_test\j01.scep';
% srcScep = loadBin(fsrcScep, 'float', 19);
% clear fsrcScep
% % f0
% fsrcF0 = 'J:\VoiceConversion\STRAIGHT-based\fws\f0_test\j01.f0';
% srcF0 = load(fsrcF0);
% clear fsrcF0
% % fout
% foutWav = 'J:\VoiceConversion\STRAIGHT-based\fws2mht_withDelta_reduced1of5\test\j01.wav';
% 
% fnum = size(Y, 2);
% Y(:, fnum) = [];
% Y(:, 1) = [];
% fnum = fnum - 2;
% 
% % synthesis sound
% energy = srcScep(1, :);
% Y      = [energy; Y(1:17, :)];
% ap     = repmat(-20, 513, fnum);
% fs     = 16000;
% sgram  = cep2sgram(Y);
% [sy, prmS] = exstraightsynth(srcF0, sgram, ap, fs);
% 
% % write to wav file
% wavwrite(sy/32768, fs, 16, foutWav);