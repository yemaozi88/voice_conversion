%
% 2011/09/28
% simple version of MLtest4MixNum.m
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%

clc, clear all, fclose('all');

%% definition
dirIn_train = 'J:\H2Swith16deg_0243\joint\jointS2H_11_1of8';
%dirIn_train = 'J:\!gesture\transitionAmong16of28\dgvs\1';
%dirIn_train = 'C:\research\H2Swith16deg_0243\!speech\ATR_A\scep16_1of8';
dirIn_test  = 'J:\H2Swith16deg_0243\joint\jointS2H_22_1of8';
%dirIn_test = 'J:\!gesture\transitionAmong16of28\dgvs\2';
%dirIn_test = 'C:\research\H2Swith16deg_0243\!speech\ATR_J\scep16_1of8';
dirOut      = 'J:\H2Swith16deg_0243\joint\H2Smodel_full_withPCA';
EigenParamDir = 'C:\research\!gesture\transitionAmong16of28\EigenParam16\1';
gMax_       = 4096;

% load Eigen parameters
[EVec, EVal, u] = loadEigenParam(EigenParamDir);
clear EigenParamDir

% % gesture
% Y = loadBinDir(dirIn_train, 'uchar', 26);
% Y = Y(5:22, :);
% Y = Y(1:16, :)';
% Y = PCA_Trans(Y, EVec, u, 16);
% Y = adddelta(Y');
% Y = Y';
% 
% U = loadBinDir(dirIn_test, 'uchar', 26);
% U = U(5:22, :);
% U = U(1:16, :)';
% U = PCA_Trans(U, EVec, u, 16);
% U = adddelta(U');
% U = U';

% % delta
% Y = Y';
% Y = adddelta(Y)';
% %Y = PCA_Trans(Y, EVec, u, 32);
% 
% U = U';
% U = adddelta(U)';
% %U = PCA_Trans(U, EVec, u, 32);

% % speech
% Y = loadBinDir(dirIn_train, 'float', 17);
% Y = Y(1:16, :);
% %Y = adddelta(Y);
% Y = Y';
%  
% U = loadBinDir(dirIn_test, 'float', 17);
% U = U(1:16, :);
% %U = adddelta(U);
% U = U';

% % joint
Y = loadBinDir(dirIn_train, 'float', 32);
U = loadBinDir(dirIn_test, 'float', 32);
% %Y = Y';
% %U = U';
% %Y = adddelta(Y)';
% %U = adddelta(U)';
  
% perform PCA
YS = Y(1:16, :);
YH = Y(17:32, :);
YHpca = PCA_Trans(YH', EVec, u, 16);
YHpca = YHpca';
% %YSsd = adddelta(YS);
% %YHsd = adddelta(YH);
% %YHsd = adddelta(YHpca);
% %clear YHpca
% % %YSd = YSsd(17:32, :);
% % %YHd = YHsd(17:32, :);
Y = [YHpca; YS];
% % Y = [YHsd; YSsd];
% % % Y = [YS; YHpca; YSd; YHd];
%  
% % the first and the last data are wrong because of delta
% % fmax = size(Y, 2);
% % Y(:, fmax) = [];
% % Y(:, 1) = [];
Y = Y';
clear YS YH YHpca
% % YSsd YHsd fmax
% 
% % % this part is used for S2HwithDelta, but it is wrong
% % % %Y = [YHpca, YS];
% % % %Y = [YS, YHpca];
% % % Y = [YSd', YHpca];
% % % %clear YH YHpca YS
% % % clear YS YH YSd YHd YHpca
 
US = U(1:16, :);
UH = U(17:32, :);
UHpca = PCA_Trans(UH', EVec, u, 16);
UHpca = UHpca';
% %USsd = adddelta(US);
% %UHsd = adddelta(UH);
% %UHsd = adddelta(UHpca);
% %clear UHpca
% % %USd = USsd(17:32, :);
% % %UHd = UHsd(17:32, :);
U = [UHpca; US];
% % U = [UHsd; USsd];
% % %U = [US; UH; USd; UHd];
% % % U = [US; UHpca; USd; UHd];
%  
% % the first and the last data are wrong because of delta
% % fmax = size(U, 2);
% % U(:, fmax) = [];
% % U(:, 1) = [];
U = U';
clear US UH UHpca
% % clear USsd UHsd fmax


%% preparation for ML calculation
gNumArray = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];
gMax = length(gNumArray);
gMax = log2(gMax_)+1;
ML = zeros(gMax, 1);
clear gMax_

% log
flog = fopen([dirOut '\log.txt'], 'wt');
fprintf(flog, '<Maximum Likelihood>\n');
fprintf(flog, 'The number of data : %d\n\n', size(Y, 1));

%for g = 1:gMax
for g = 4:4   
tic
    gNum = gNumArray(g);
    disp(gNum)

    obj = trainGMM(Y, gNum, 0); % 0 - full, 1 - diagonal
    fOut = [dirOut '\H2Smodel_mix' num2str(gNum) '_obj'];
    save('obj');
    %load(fOut)
    clear fOut
        
    L_ = pdf(obj, U);
    
%     % remove Inf and -Inf frames
%     fmax_ = length(L_);
%     L = [];
%     for ii = 1:fmax_
%         if log(L_(ii, 1)) ~= Inf && log(L_(ii, 1)) ~= -Inf
%             L = [L; log(L_(ii, 1))];
%         end
%     end
%     fmax = length(L);
%     cutoff = fmax_ - fmax;
%     disp([num2str(cutoff) ' samples are removed.'])
%     L_ = L;
%     clear L fmax fmax_
%     
%     L  = mean(L_);

    L = log(max(L_));
    ML(g, 1) = L;
    disp([num2str(gNum) ':' num2str(L)]);

    % log
    fprintf(flog, 'The number of mixtures : %d\n', gNum);
    %fprintf(flog, '\ttrain : %6.4f\n', obj.NlogL);
    fprintf(flog, '\ttest  : %6.4f\n', L);
    
    save([dirOut '\ML'], 'ML');
toc
end
fclose(flog);
clear flog