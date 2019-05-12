function Y = gmmMappingWithDelta2(X, obj)
% Y = gmmMappingWithDelta(X, obj)
% estimates Y using obj trained by parallel data of [X, Y] which include
% delta features (static part and delta part are treated at the same time)
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
% %fin = 'J:\VoiceConversion\STRAIGHT-based\fws\scep18+delta_test\j01.feature';
% %X = loadBin(fin, 'float', 36);
% % hand motion
% fin = 'J:\H2SwithDelta\!dgvd\1\ai.dgvd';
% X = loadBin(fin, 'float', 36);
% clear fin
%                      
% % GMM object
% %load('C:\research\VoiceConversion\STRAIGHT-based\fws2mht_withDelta_reduced1of5\gmm_full\jointModel_mix8_obj');
% load('J:\H2SwithDelta\withDelta\gmm_full\jointModel_mix4_obj');


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
muY = obj.mu(:, y_idx); % gNum x tDim, (tDim = jDim - sDim)
p = obj.PComponents; % NComponents x 1, proportions

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
end
%clear x_idx y_idx obj

    

%% E and D matrice
% Y = zeros(size(X));
objX = gmdistribution(muX, sigmaXX, p);

for ii = 1:sNum
    x = X(:, ii);

    % choose the most likely distribution
    pX = posterior(objX, x'); % X should be fNum x dim
    [maxValue, jj] = max(pX);
    clear pX maxValue
    
    % sDim x 1
    sigmaXX_ = squeeze(sigmaXX(:, :, jj))';
    sigmaXY_ = squeeze(sigmaXY(:, :, jj))';
    sigmaYX_ = squeeze(sigmaYX(:, :, jj))';
    sigmaYY_ = squeeze(sigmaYY(:, :, jj))';

    E(:, ii) = muY(jj, :)' + sigmaYX_ ./ sigmaXX_ .* (x - muX(jj, :)'); % sDim x 1
    D(:, ii) = sigmaYY_ - sigmaYX_ ./ sigmaXX_ .* sigmaXY_; % sDim x 1
    
    clear sigmaXX_ invSigmaXX_ sigmaXY_ sigmaYX_ sigmaYY_
end %sNum
clear x ii jj s_idx d_idx
clear muX muY p objX
clear sigmaXX sigmaXY sigmaYX sigmaYY

Evec = reshape(E, sDim * sNum, 1);
Dvec = reshape(D, sDim * sNum, 1);
clear E D


%% W matrix
row = sNum * sDim;
col = sNum * dDim;

% weight: static: 1, delta: 0.5
W = spalloc(row, col, 3*col);
for t = 1:sNum
    for h = 1:dDim
        % static
        W(2*(t-1)*dDim + h, (t-1)*dDim + h) = 1;
        % delta
        if t == 1
            W((2*t-1)*dDim + h, t*dDim + h) = 0.5;
        elseif t == sNum
            W((2*t-1)*dDim + h, (t-2)*dDim + h) = -0.5;
        else
            W((2*t-1)*dDim + h, (t-2)*dDim + h) = -0.5;
            W((2*t-1)*dDim + h, t*dDim + h) = 0.5;
        end
    end % h 
end % t
clear t h row col


%% mapping
Dm = spalloc(sDim*sNum, sDim*sNum, sDim*sNum);
for mDim = 1 : sDim*sNum
     Dm(mDim, mDim) = Dvec(mDim);
end
invDm = inv(Dm);
clear mDim Dvec Dm;

% eq. 39
traW_invDm = W' * invDm;
A = traW_invDm * W;
b = traW_invDm * Evec;
clear invDm traW_invDm W Evec

Y = A \ b;
Y = reshape(Y, dDim, length(Y)/dDim);
clear A b
clear sDim jDim gNum sNum


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
% %Y(:, fnum) = [];
% %Y(:, 1) = [];
% %fnum = fnum - 2;
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