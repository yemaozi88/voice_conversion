%
% 2011/07/09
% S2H-H2S combined system to choose best model for F0 generation
% gesture data is performed PCA
% speech data is added F0 and VS, then performed PCA
%
% LINK
% makeGestureCombinationList.m
%
% HISTORY 
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%

clc, clear all, fclose('all');


%% definition
% research folder
drive = 'D:'; % it can be 'C:\research' or 'G:'

% gesture/speaker model
fGestureModel = [drive '\F0generation\objGestureModel-32_withPCA'];
fSpeakerModel = [drive '\F0generation\preliminaryTest\1113_withPCA\objSpeakerModel-128'];

% PCA info
EigenParamDirH = [drive '\!gesture\transitionAmong16of28\EigenParam\1'];
EigenParamDirS = [drive '\F0generation\!speech\ATR_A\EigenParam_ATR_A'];

% gesture and speech for convert model
dirH = [drive '\!gesture\transitionAmong16of28\dgvs'];
dirS = [drive '\F0generation\!speech\vowels\feature_scep-f0-vs'];

% the directory for the output files
dirOut = [drive '\F0generation\8190combinations_withPCA'];

% S2H model, full covariance
MIX = 8;

% gesture and speech for convertion
dirVowels      = [drive '\F0generation\!speech\vowels\feature_scep-f0-vs\1'];
dirConsonants_ = [drive '\F0generation\!speech'];

consonant = ['b', 'm', 'n', 'p', 'r'];
vowel = ['a', 'i', 'u', 'e', 'o'];

nS = 1; % the number of speech dataset

% Saito's method
alpha = 1; % weight factor for speaker model
it = 1; % number or iteration
updatemethod = 1; %0- using target responsibility 1- using joint responsibility
 
SAMPLING_FREQ = 1; % assumed sampling frequency of DataGlove

% feature to wav conversion
WINDOW = 40;  % window size for smoothing
THRES  = 0.75; % threshold for vocing strength

% cepstralDistance
%cepstralDistancefile = [dirOut '\cepstralDistortion.csv'];
clear drive


%% Gesture model / Speaker model
% covariance matrice are diagonal
load(fGestureModel);
load(fSpeakerModel);        
clear fGestureModel fSpeakerModel

% load PCA info
[EVecH, EValH, EuH] = loadEigenParam(EigenParamDirH);
[EVecS, EValS, EuS] = loadEigenParam(EigenParamDirS);
clear EigenParamDirH EigenParamDirS


%% choose 5 gestures among 16 gestures
ges = [1, 2, 4, 7, 8, 9, 11, 13, 14, 15, 16, 21, 22, 25, 27, 28];
[dgvMean, gestureCombinationList] = make16gestureCombinationList(dirH);
P = gestureCombinationList;
Pnum = size(P, 1);  % the number of all permutations

%fcepDis = fopen(cepstralDistancefile, 'wt');

for nP = 8000:8192;  % the number of the permutation
tic
    disp('--- make joint vectors ---')
    tic
    % set vowel gestures
    
    a = P(nP, 1);
    if a < 10
        a = num2str(a);
        a = ['0' a];
    else
        a = num2str(a);
    end

    i = P(nP, 2);
    if i < 10
        i = num2str(i);
        i = ['0' i];
    else
        i = num2str(i);
    end

    u = P(nP, 3);
    if u < 10
        u = num2str(u);
        u = ['0' u];
    else
        u = num2str(u);
    end

    e = P(nP, 4);
    if e < 10
        e = num2str(e);
        e = ['0' e];
    else
        e = num2str(e);
    end

    o = P(nP, 5);
    if o < 10
        o = num2str(o);
        o = ['0' o];
    else
        o = num2str(o);
    end

    %% nP: num 2 str
    if nP < 10
        nPstr = ['000' num2str(nP)];
    elseif nP < 100
        nPstr = ['00' num2str(nP)];
    elseif nP < 1000
        nPstr = ['0' num2str(nP)];
    else
        nPstr = num2str(nP);
    end
    disp(nPstr)

    
%% directory processing
    if ismac == 1
        dirOutSub          = [dirOut '/' num2str(nPstr)];
        dirSynDgv          = [dirOutSub '/synDgv'];
        dirSynFeature      = [dirOutSub '/synFeature'];
%        dirOutJoint           = [dirOutSub '/joint'];
%        dirOutReducedJoint    = [dirOutSub '/reducedJoint'];
%        dirOutReducedJointLog = [dirOutReducedJoint '\log'];
    else
        dirOutSub          = [dirOut '\' num2str(nPstr)];
        dirSynDgv          = [dirOutSub '\synDgv'];
        dirSynFeature      = [dirOutSub '\synFeature'];
%        dirOutJoint           = [dirOutSub '\joint'];
%        dirOutReducedJoint    = [dirOutSub '\reducedJoint'];
%        dirOutReducedJointLog = [dirOutReducedJoint '\log'];
    end    
    mkdir(dirSynDgv);
    mkdir(dirSynFeature);
%    mkdir(dirOutJoint);
%    mkdir(dirOutReducedJoint);
%    mkdir(dirOutReducedJointLog);
    

%% log
if ismac == 1
    fname_log = [dirOutSub '/log.txt'];
else
    fname_log = [dirOutSub '\log.txt'];
end
flog  = fopen(fname_log, 'wt');
fprintf(flog, '< Gestures for the five Japanese vowels >\n');

fprintf(flog, 'a: %s\t', a);
fprintf(flog, 'i: %s\t', i);
fprintf(flog, 'u: %s\t', u);
fprintf(flog, 'e: %s\t', e);
fprintf(flog, 'o: %s\n\n', o);

disp(['Gesture Design: ' num2str(nP)]);
disp(['a:' a ' i:' i ' u:' u ' e:' e ' o:' o]);


    %% load files
    Jnum = 0;
    Y = [];
    for nH = 1:1; % the number of dgvs data set
        nS = nH;
% log
fprintf(flog, '< The number of frames >\n');
fprintf(flog, 'dataset: Hand %d, Speech %d\n', nH, nS);

        for nV1 = 1:5 % the first vowel of a transition in training data
            for nV2 = 1:5 % the second vowel of a transition in training data

% log
fprintf(flog, '%s%s:\t', vowel(nV1), vowel(nV2));

                % load gesture data
                if ismac == 1
                    filenameH_ = sprintf('filenameH = [''%s/'' ''%d'' ''/'' %s ''-'' %s ''.dgvs''];', dirH, nH, vowel(nV1), vowel(nV2));
                else
                    filenameH_ = sprintf('filenameH = [''%s\\'' ''%d'' ''\\'' %s ''-'' %s ''.dgvs''];', dirH, nH, vowel(nV1), vowel(nV2));
                end

                eval(filenameH_);
                H = loadBin(filenameH, 'uchar', 26);
                frameH = size(H, 2);
                clear filenameH_

                % load speech feature data
                if ismac == 1
                    filenameS = [dirS '/' num2str(nS) '/' vowel(nV1) vowel(nV2) '.feature'];
                else
                    filenameS = [dirS '\' num2str(nS) '\' vowel(nV1) vowel(nV2) '.feature'];
                end
                S = loadBin(filenameS, 'float', 18);
                frameS = size(S, 2);

% log
fprintf(flog, 'gesture,%4d\t', frameH);
fprintf(flog, 'speech,%4d\t', frameS);

                % make augmented vector
                J = makeJoint(H, S, 1, 0); % with scep0, J = [H, S]
                frameJ = size(J, 2);
                Jnum = Jnum + frameJ;

                Y = [Y, J];
                
                % save augmented vector
%                 if ismac == 1
%                     fname_joint = [dirOutReducedJointSub '/' vowel(nV1) vowel(nV2) '_' num2str(nH) num2str(nS) '.joint'];
%                 else
%                     fname_joint = [dirOutReducedJointSub '\' vowel(nV1) vowel(nV2) '_' num2str(nH) num2str(nS) '.joint'];
%                 end
%                 fjoint = fopen(fname_joint, 'wb');
%                 for ii = 1:frameK
%                     fwrite(fjoint, K(:, ii), 'float');
%                 end
%                 fclose(fjoint);

% log
fprintf(flog, 'joint,%4d\n', frameJ);

                %clear fname_joint fjoint frameJ ii;
                clear J frameJ;
            end % nV2
        end % nV1
% log
fprintf(flog, '\n');

        %end % nS
    end % nH
%log
fprintf(flog, 'total frame number in all joint vectors: %4d\n', Jnum);
    fclose(flog);
    clear flog fname_log
    clear filenameH filenameS
    clear frameH frameS
    clear nV1 nV2 nH
    clear H S J Jnum   


%% GMM training
    disp('--- GMM training ---')
    Y = Y';
    
    % for HS joint data
    YH = Y(:, 1:18);
    YS = Y(:, 19:36);
    YH = PCA_Trans(YH, EVecH, EuH, 18);
    YS = PCA_Trans(YS, EVecS, EuS, 18);
    Y  = [YH, YS];
    clear YH YS
    
    objH2Smodel = trainGMM(Y, MIX, 0);
    clear Y
    
    if ismac == 1
        save([dirOutSub '/objH2Smodel'], 'objH2Smodel');
    else
        save([dirOutSub '\objH2Smodel'], 'objH2Smodel');
    end


%% S2H
    disp('--- S2H conversion ---')
    
    %CD = zeros(30, nSmax); % Cepstral Distortion for vowels + consonants
    for nC = 1:5; % consonant
        for nV1 = 0:5 % 0 means consonant
            if nC == 1 || nV1 == 0 

                for nV2 = 1:5
                    if nV1 == 0             
                        mora  = sprintf('%s%s', consonant(nC), vowel(nV2));
                        fname = mora;
                        dirConsonants = [dirConsonants_ '\' consonant(nC) '\feature_scep-f0-vs\' num2str(nS)];

                        if ismac == 1
                            fname_scepIn = [dirConsonants '/' fname '.feature'];
                            fname_dgvSyn = [dirSynDgv '/' fname '.dgv'];
                            fname_dgvLog = [dirSynDgv '/' fname '.txt'];
                        else
                            fname_scepIn = [dirConsonants '\' fname '.feature'];
                            fname_dgvSyn = [dirSynDgv '\' fname '.dgv'];
                            fname_dgvLog = [dirSynDgv '\' fname '.txt'];
                        end

                    else
                        mora  = sprintf('%s%s', vowel(nV1), vowel(nV2));
                        fname = mora;
                        if ismac == 1
                            fname_scepIn = [dirVowels '/' fname '.feature'];
                            fname_dgvSyn = [dirSynDgv '/' fname '.dgv'];
                            fname_dgvLog = [dirSynDgv '/' fname '.txt'];
                        else
                            fname_scepIn = [dirVowels '\' fname '.feature'];
                            fname_dgvSyn = [dirSynDgv '\' fname '.dgv'];
                            fname_dgvLog = [dirSynDgv '\' fname '.txt'];
                        end
                    end
                    %disp(mora)
        
                    % conversion
                    input = loadBin(fname_scepIn, 'float', 18);
                    input = input';
                    input = PCA_Trans(input, EVecS, EuS);
                    input = input';
        
                    dgv_ = spkmodel_vc2_(input, objH2Smodel, objGestureModelWithPCA, alpha, it, updatemethod, fname_dgvLog);
                    %dgv_ = spkmodel_vc2(input, objS2Hmodel, objGestureModelWithPCA, alpha, it, updatemethod, fname_dgvLog);
                    %dgv_ = gmmvc(input, objS2Hmodel);
                    % dgv_ holds all results at every step
                    dgv_ = dgv_{1, it};
        
                    % perform invPCA
                    dgv_ = dgv_';
                    dgv_ = PCA_TransInv(dgv_, EVecH, EuH);
                    dgv_ = dgv_';
        
                    dgv  = conv2dgv(dgv_, SAMPLING_FREQ);
                    frameDgv = size(dgv, 2);
        
                    % write out
                    fout = fopen(fname_dgvSyn, 'wb');
                    for ii = 1:frameDgv
                        fwrite(fout, dgv(:, ii), 'uchar');
                    end
                    fclose(fout);
                    clear input dgv dgv_ frameDgv
                end % nV2
            end % if
        end % nV1
    end % nC
    clear fname_dgvLog fname_dgvSyn fname_scepIn
    
    
%% H2S
    disp('--- H2S conversion ---')
    
    dirlist = dir([dirSynDgv '\*.dgv']);
    dirlength = length(dirlist);

    ii = 1;
    while ii < dirlength + 1
        % except ".", "..", "DS_Store"
        if length(dirlist(ii).name) > 3 
            filename = strrep(dirlist(ii).name, '.dgv', '');
            %disp(filename)

            if ismac == 1
                fname_dgvs      = [dirSynDgv '/' filename '.dgv'];
                fname_scep      = [dirSynFeature '/' filename '.feature'];
                fname_scepLog   = [dirSynFeature '/' filename '.txt'];
            else
                fname_dgvs      = [dirSynDgv '\' filename '.dgv'];
                fname_scep      = [dirSynFeature '\' filename '.feature'];
                fname_scepLog   = [dirSynFeature '\' filename '.txt'];
            end

            input = loadBin(fname_dgvs, 'uchar', 26);
            input = input(5:22, :);

            % perform PCA
            input = input';
            input = PCA_Trans(input, EVecH, EuH);
            input = input';
            
            scep_ = spkmodel_vc2(input, objH2Smodel, objSpeakerModel, alpha, it, updatemethod, fname_scepLog);
            %scep_ = spkmodel_vc2_(input, objS2Hmodel, objSpeakerModel, alpha, it, updatemethod, fname_scepLog);
            %scep_ = gmmvc(input, objH2Smodel);
            % scep_ holds all results at every step
            scep_ = scep_{1, it};
            
            % perform invPCA
            scep_ = scep_';
            scep_ = PCA_TransInv(scep_, EVecS, EuS);
            scep  = scep_';
            clear scep_

            % write out
            frameScep = size(scep, 2);
            fod = fopen(fname_scep, 'wb');
            j = 1;
            while j < frameScep + 1
                fwrite(fod, scep(:, j), 'float');
                j = j + 1;
            end
            fclose(fod);
            clear scep
            
            % feature to wav
            feature2wav(fname_scep, WINDOW, THRES);
        end
        ii = ii + 1;
    end % while
    toc
    clear ii j nC nV1 nV2
    clear fname fname_dgvs fname_scep fname_scepLog
    

%% comparison
%     %CDmean = mean(CD')'; % 30 x 1
%     CDmean = CD;
%     CDmax = size(CDmean, 1);
%     
%     fprintf(fcepDis, '%d,%s,%s,%s,%s,%s,', nP, a, i, u, e, o);
%     for ii = 1:CDmax
%         fprintf(fcepDis, '%f,', CDmean(ii, 1));
%     end
%     fprintf(fcepDis, '\n');
% toc
end % nP