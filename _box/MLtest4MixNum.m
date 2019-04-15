%
% 2010/11/25
% MLtest4MixNum.m evaluates the appropriate mixture number calculating ML
% 
% NOTES
% - followings are unchecked
%   TYPE: gesture, speech
%   TEST: 1, 3
%
% HISTORY
% 2011/06/21 united 3 parts, gesture, speech and joint part
%
% AUTHOR
% Aki Kunikoshi (D2)
% yemaozi88@gmail.com
%

clc, clear all, fclose('all');

%% definition
dirIn_train = 'K:\!gesture\transitionAmong16of28\dgv\1';
dirIn_test  = 'K:\!gesture\transitionAmong16of28\dgv\2';
dirOut      = 'K:\realtimeDemo\v6\demoSamples\gestureModel_full';

%EigenParamDirH = 'J:\!gesture\transitionAmong16of28\EigenParam\1';
%EigenParamDirS = 'J:\F0generation\withPCA\EigenParam_ATR_A';
%EigenParamDirH = 'C:\research\H2SwithDelta\!dgvd\EigenParam\1';
%EigenParamDirS = 'C:\research\H2SwithDelta\!scepd\EigenParam\1';

% load EigenParam
%[EVecH, EValH, EuH] = loadEigenParam(EigenParamDirH);
%[EVecS, EValS, EuS] = loadEigenParam(EigenParamDirS);
%clear EigenParamDirH

% TYPE
%   0: gesture (26 uchar)
%   1: speech (18 float)
%   2: joint (36 float)
%   3: feature (scep 1-16 + F0 + VS, 18 float)
%   4: feature joint (72 float)
%   5: speech + delta (36 float)
%   6: gesture + delta (36 float)
TYPE = 0;

% TRAIN
%   0: gmm train and save obj
%   1: load obj
TRAIN = 0;

% TEST
%   0: each file in the directory
%   1: specific files in the directory
%   2: all files in the directory
%   3: else
TEST = 2;
Smax = 25;

% Maximum mixture number
gMax_ = 1024;

% output file name
logfile = 'log_all.txt';
MLfile  = 'ML_all';


%% load training data
if TYPE == 0; % gesture
    disp('data type: gesture')
    Y = loadBinDir(dirIn_train, 'uchar', 26);
elseif TYPE == 1; % speech
    disp('data type: cepstrum')
    Y = loadBinDir(dirIn_train, 'float', 18);
elseif TYPE == 2; % joint
    disp('data type: joint')
    Y = loadBinDir(dirIn_train, 'float', 36);
elseif TYPE == 3; % feature (scep 1-16 + F0 + VS, 18 float)
    disp('data type: feature (scep 1-16 + F0 + VS)')
    Y = loadBinDir(dirIn_train, 'float', 18);
elseif TYPE == 4; % feature joint (72 float)
    disp('data type: feature joint')
    Y = loadBinDir(dirIn_train, 'float', 72);
elseif TYPE == 5 || TYPE == 6; % scep + delta (36 float)
    disp('data type: dgv/scep + delta')
    Y = loadBinDir(dirIn_train, 'float', 36);    
end

%YH = Y(1:18, :)';
%YS = Y(37:54, :)';

% perform PCA

% YH = Y(:, 1:18);
% YS = Y(:, 19:36);
% YH = PCA_Trans(YH, EVecH, EuH, 18);
% YS = PCA_Trans(YS, EVecS, EuS, 18);
% Y  = [YH, YS];
% clear YH YS

% Ydgv  = Y(1:36, :)';
% Yscep = Y(37:72, :)';
% YdgvPCA  = PCA_Trans(Ydgv, EVecH, EuH, 36);
% YscepPCA = PCA_Trans(Yscep, EVecS, EuS, 36);

%Y = [YdgvPCA, Yscep];
%Y(:, 70:72) = Y(:, 70:72) * 10; % covariance got too small...
% clear Ydgv Yscep YdgvPCA

%Y = Y';

% gesture
Y = Y(5:20, :)';

%YHpca = PCA_Trans(YH, EVecH, EuH, 18);
%Y = [YHpca, YS];
%clear YH YHpca YS


%% preparation for ML calculation
gNumArray = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
gMax = length(gNumArray);
gMax = log2(gMax_)+1;
clear gMax_

% log
flog = fopen([dirOut '\' logfile], 'wt');
fprintf(flog, '<Maximum Likelihood>\n');
fprintf(flog, 'The number of data : %d\n\n', size(Y, 1));

if TRAIN == 1;
    clear Y
end

if TEST == 0; % each file in the directory
%    gMax = 9;
    ML = zeros(gMax, Smax);

    % index of each data
%     idx = ['aa'; 'ae'; 'ai'; 'ao'; 'au'; ...
%         'ea'; 'ee'; 'ei'; 'eo'; 'eu'; ...
%         'ia'; 'ie'; 'ii'; 'io'; 'iu'; ...
%         'oa'; 'oe'; 'oi'; 'oo'; 'ou'; ...
%         'ua'; 'ue'; 'ui'; 'uo'; 'uu'; 'MN'];
%     %idx = cellstr(idx); % 25x1 cell

% log
% for ii = 1:Smax+1;
%     fprintf(flog, '%s\t',idx(ii, :));
% end
%     clear idx
%     fprintf(flog, '\n');

elseif TEST == 1; % specific files in the directory
    ML = zeros(gMax, 1);
elseif TEST == 2; % all files in the directory
    ML = zeros(gMax, 1);
else
end


%% calculate ML
for g = 9:gMax
    tic
    gNum = gNumArray(g);
    disp(gNum)

    % get object
    if TRAIN == 0;
        obj = trainGMM(Y, gNum, 0); % 0 - full, 1 - diagonal
        %[xCovArr, xMeanArr, WArr, gNum] = ConvertMatlabGMM(obj);
        if TYPE == 1 || TYPE == 5;
            fOut = [dirOut '\SpkModel_mix' num2str(gNum) '_obj'];
        elseif TYPE == 2 || TYPE == 4;
            fOut = [dirOut '\jointModel_mix' num2str(gNum) '_obj'];
        elseif TYPE == 0 || TYPE == 6;
            fOut = [dirOut '\gestureModel_mix' num2str(gNum) '_obj'];           
        end

        save(fOut, 'obj');
        clear fOut
    elseif TRAIN == 1;
        if TYPE == 1 || TYPE == 5;
            load([dirOut '\SpkModel_mix' num2str(gNum) '_obj']);
        elseif TYPE == 2 || type == 4;
            load([dirOut '\jointModel_mix' num2str(gNum) 'obj']);
        elseif TYPE == 0 || TYPE == 6;
            load([dirOut '\gestureModel_mix' num2str(gNum) '_obj']);        
        end
    end
 
    % ML calculation

    % TEST
    %   0: each file in the directory
    %   1: specific files in the directory
    %   2: all files in the directory
    %   3: else
    
%%%%%% TEST == 0 %%%%%
    if TEST == 0; % each file in the directory
        dirlist = dir(dirIn_test);
        dirlength = length(dirlist);

        Snum = 1; % the number of samples in the test set
        for ii = 1:dirlength
            % except ".", "..", "DS_Store"
            if length(dirlist(ii).name) > 3
                filename = dirlist(ii).name;
                if ismac == 1
                    fin = [dirIn_test '/' filename];
                else
                    fin = [dirIn_test '\' filename];
                    if TYPE == 0; % gesture
                        U = loadBin([dirIn_test '\' filename], 'uchar', 26);
                    elseif TYPE == 1; % speech
                        U = loadBin([dirIn_test '\' filename], 'float', 19);
                    elseif TYPE == 2; % joint
                        U = loadBin([dirIn_test '\' filename], 'float', 36);
                    elseif TYPE == 3; % feature (scep 1-16 + F0 + VS, 18 float)
                        U = loadBin([dirIn_test '\' filename], 'float', 18);
                    elseif TYPE == 4; % feature joint
                        U = loadBin([dirIn_test '\' filename], 'float', 72);
                    elseif TYPE == 5; % scep + delta
                        U = loadBin([dirIn_test '\' filename], 'float', 36);
                    end
                end

                %U = U';

                % perform PCA

%                 UH = U(:, 1:18);
%                 US = U(:, 19:36);
%                 UH = PCA_Trans(UH, EVecH, EuH, 18);
%                 US = PCA_Trans(US, EVecS, EuS, 18);
%                 U  = [UH, US];
%                 clear UH US

                %Udgv  = U(1:36, :)';
                %Uscep = U(37:72, :)';
 
                %UdgvPCA  = PCA_Trans(Udgv, EVecH, EuH, 36);
%                 UscepPCA = PCA_Trans(Uscep, EVecS, EuS, 36);

                %U = [UdgvPCA, Uscep];
%                 U(:, 70:72) = U(:, 70:72) * 10; % covariance got too small...
                %clear Udgv Uscep UdgvPCA

                UH = U(1:18, :)';
                US = U(37:54, :)';
                UHpca = PCA_Trans(UH, EVecH, EuH, 18);
                U = [UHpca, US];
                clear UH UHpca US


                L_ = pdf(obj, U);
                L  = log(max(L_)); % maximum likelihood
                disp([num2str(gNum) ' - ' filename ':' num2str(L)]);
% log
%fprintf(flog, 'test-%s, %6.4f\n',filename, L);
                ML(g, Snum) = L;
                Snum = Snum + 1;
            end
        end
        %clear dirlist dirlength filename
        %clear ii fin L L_ U

%%%%%% TEST == 1 %%%%%
    elseif TEST == 1
        % specific files in the directory
        MLavg = 0;
        nnMax = 10;
        for nn = 1:nnMax;
            if nn < 10
                fnameTest = [dirIn_test '\j0' num2str(nn) '.joint'];
            else
                fnameTest = [dirIn_test '\j' num2str(nn) '.joint'];
            end
            U = loadBin(fnameTest, 'float', 36);
            U = U';
            %U = U(2:19, :)';
            %U = U(1:18, :)';
            
            % perform PCA
%             UH = U(:, 1:18);
%             US = U(:, 19:36);
%             UH = PCA_Trans(UH, EVecH, EuH, 18);
%             US = PCA_Trans(US, EVecS, EuS, 18);
%             U  = [UH, US];
%             clear UH US
            
            L_ = pdf(obj, U);
            L  = log(max(L_)); % maximum likelihood            

% log
fprintf(flog, 'The number of mixtures : %d\n', gNum);
fprintf(flog, 'train, %6.4f\n', obj.NlogL);
if nn < 10
    fprintf(flog, 'test-0%d, %6.4f\n', nn, L);
else
    fprintf(flog, 'test-%d, %6.4f\n', nn, L);
end

            MLavg = MLavg + L;
            clear fnameTest L_ L U
        end % nn
        MLavg = MLavg / nnMax;
        disp(['Average : ' num2str(MLavg)])
% log
fprintf(flog, 'Average, %6.4f\n\n', MLavg);
            
        ML(g, 1) = MLavg;

%%%%%% TEST == 2 %%%%%
    elseif TEST == 2
        % all files in the directory
        if ismac == 1
        else
            if TYPE == 0; % gesture
                U = loadBinDir(dirIn_test, 'uchar', 26);
            end
        end
        %U = U';

%         % perform PCA
%         UH = U(:, 1:18);
%         US = U(:, 19:36);
%         UH = PCA_Trans(UH, EVecH, EuH, 18);
%         US = PCA_Trans(US, EVecS, EuS, 18);
%         U  = [UH, US];
%         clear UH US

%         U = U';
%         Upca = PCA_Trans(U, EVecH, EuH, 36);
%         U = Upca;
%         clear Upca
        U = U(5:20, :)';

        L_ = pdf(obj, U);
        L  = log(max(L_)); % maximum likelihood
        disp([num2str(gNum) ' : ' num2str(L)]);
% log
fprintf(flog, 'The number of mixtures : %d\n', gNum);
fprintf(flog, 'train, %6.4f\n', obj.NlogL);
fprintf(flog, 'Average, %6.4f\n\n', L);        
        ML(g, 1) = L;
        %clear U L_ L
    else
            % %     %DirIn_ = 'L:\research\gesture\randomMove20\dgvs\';
% %     %ML_ = zeros(10, 1);
% % 
% % %     %% gesture data
% % %     % for training
% % %     for j = 0:9
% % %         X = [];
% % %         j_ = num2str(j);
% % %         DirIn = [DirIn_ 'rand' j_ '.dgvs'];
% % %         %disp(DirIn);
% % %         T = loadBin(DirIn, 'uchar', 26);
% % %
% % %         % for testing (cross validation)
% % %         for i = 0:9
% % %             if i ~= j
% % %                 i_ = num2str(i);
% % %                 DirIn = [DirIn_ 'rand' i_ '.dgvs'];
% % %                 %disp(DirIn);
% % %                 tmp = loadBin(DirIn, 'uchar', 26);
% % %                 X = [X, tmp];
% % %             end
% % %         end
    end
             
    save([dirOut '\' MLfile], 'ML');
    toc
end % g
%clear g gNum obj U


%% add mean
if TEST == 0

% %% non sence to calculate mean here!
%     MLmeanG = mean(ML')'; % mean over gNum
%     ML = [ML, MLmeanG];
% 
%     MLmeanS = sum(ML, 1); % mean for every sample
%     MLmeanS = MLmeanS ./ gMax;
%     ML = [ML; MLmeanS];
% 
%     clear MLmeanG MLmeanS

    % log
    for g = 1:gMax
        gNum = gNumArray(g);
        for Snum = 1:Smax;
            fprintf(flog, '%6.2f\t', ML(g, Snum));
        end
        fprintf(flog, '\n');
    end
    %clear g gNum gNumArray
    %clear Smax Snum
end

fclose(flog);
clear flog


%% gesture model
% training data
%X = loadBinDir('C:\research\_gesture\randomMove\dgvs', 'uchar', 26);
%Y = loadBinDir('C:\research\_gesture\transitionAmong16of28\dgvs\1', 'uchar', 26);

% evaluation data
%U1 = loadBinDir('C:\research\_gesture\transitionAmong16of28\dgvs\2', 'uchar', 26);
%U2 = loadBinDir('C:\research\_gesture\transitionAmong16of28\dgvs\3', 'uchar', 26);

%Y = Y(5:22, :)';
%U1 = U1(5:22, :)';
%U2 = U2(5:22, :)';

% eigen parameters
% X1 = loadBinDir('C:\research\_gesture\transitionAmong16of28\dgvs\1', 'uchar', 26);
% X2 = loadBinDir('C:\research\_gesture\transitionAmong16of28\dgvs\2', 'uchar', 26);
% X3 = loadBinDir('C:\research\_gesture\transitionAmong16of28\dgvs\3', 'uchar', 26);
% X = [X1, X2, X3];
% clear X1 X2 X3
% X = X(5:22, :);
% X = X';
% dirOut = 'C:\research\_gesture\transitionAmong16of28\EigenParam\all';
% getEigenParam(X, dirOut);

% EigenParamDir = 'C:\research\_gesture\transitionAmong16of28\EigenParam\1';
% [EVec, EVal, Eu] = loadEigenParam(EigenParamDir);
% Y = PCA_Trans(Y, EVec, Eu);
% U1 = PCA_Trans(U1, EVec, Eu);
% U2 = PCA_Trans(U2, EVec, Eu);


%% joint vector
%dirOutJoint_train = 'I:\VoiceConversion\mht2fws\joint_train';
%dirOutJoint_train  = 'j:\ProbabilisticIntegrationModel\distortionWithScep0\1113\0-17\reducedJoint11_HS';
%dirOutJoint_test  = 'I:\VoiceConversion\mht2fws\joint_test';
%dirOutJoint_test  = 'j:\ProbabilisticIntegrationModel\distortionWithScep0\1113\0-17\joint22_HS';

% training data
%Y = loadBinDir(dirOutJoint_train, 'float', 36);

%Y_H1S1 = loadBinDir([dirOutJoint '\H1S1'], 'float', 36);
%Y_H1S2 = loadBinDir([dirOutReducedJoint '\H1S2'], 'float', 36);
%Y_H1S3 = loadBinDir([dirOutReducedJoint '\H1S3'], 'float', 36);
%Y_H2S1 = loadBinDir([dirOutReducedJoint '\H2S1'], 'float', 36);
%Y_H2S2 = loadBinDir([dirOutReducedJoint '\H2S2'], 'float', 36);
%Y_H2S3 = loadBinDir([dirOutReducedJoint '\H2S3'], 'float', 36);

% evaluation data
%U = loadBinDir(dirOutJoint_test, 'float', 36);
%U = U';
%Y_H3S3 = loadBinDir([dirOutReducedJoint '\H3S2'], 'float', 36);

%Y = [Y_H1S1]';
%Y = [Y_H1S1, Y_H1S2]'; % 2a
%Y = [Y_H1S1, Y_H2S2]'; % 2b
%Y = [Y_H1S1, Y_H1S2, Y_H2S1, Y_H2S2]'; 


%U = Y_H3S3';
%clear Y_H1S1 Y_H1S2 Y_H2S1 Y_H2S2 Y_H3S3

% training data
% train4demo;
% dirC = 'I:\ProbabilisticIntegrationModel\distortion\5840\joint\consonant\1';
% 
% % evaluation data
% Ua = loadBin([dirC '\na-v.joint'], 'float', 36);
% Ui = loadBin([dirC '\ni-v.joint'], 'float', 36);
% Uu = loadBin([dirC '\nu-v.joint'], 'float', 36);
% Ue = loadBin([dirC '\ne-v.joint'], 'float', 36);
% Uo = loadBin([dirC '\no-v.joint'], 'float', 36);
% Ua = Ua';
% Ui = Ui';
% Uu = Uu';
% Ue = Ue';
% Uo = Uo';


%% speaker model
% % training set
% dirScep = 'F:\ProbabilisticIntegrationModel\distortion\scep18_withoutSil_a';
% 
% % log
% flog = fopen([dirOut '\dataNum4reducedScep.txt'], 'wt');
% fprintf(flog, '<The number of data>\n');
% 
% %Y = loadBinDir(dirScep, 'float', 19);
% Y = [];
% for nn = 1:50
%     if nn < 10
%         nnStr = ['0' num2str(nn)];
%     else
%         nnStr = num2str(nn);        
%     end
%     fnameScep = [dirScep '\a' nnStr '.scep'];
%     Y_ = loadBin(fnameScep, 'float', 19);
%     fnumBefore = size(Y_, 2);
%     %[Y_, fnumAfter, ratio] = reduceScep(Y_, 0.015);
%     [Y_, fnumAfter, ratio] = reduceScep(Y_, 0.018);
% 
% % log
% fprintf(flog, '%s : %4d to %4d (%4.2f[%%])\n', nnStr, fnumBefore, fnumAfter, ratio);
% disp([nnStr ' : ' num2str(fnumBefore) ' to ' num2str(fnumAfter) ' (' num2str(ratio) '[%])' ])
% 
%     Y = [Y, Y_];
% end
% clear nnStr fnameScep Y_ fnumBefore fnumAfter ratio
% fclose(flog);
% %Y = Y(2:19, :)';
% Y = Y(1:18, :)';