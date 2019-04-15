%
% 2011/07/05
% f0generation4test.m trains GMM, S2H and H2S converters with f0 generation
% speech feature is scep 1-16 deg + f0 + vs
%
% REMARKS
% this program is based on combinedModelWithScep0.m
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%

clear all, clc, fclose('all');

mode = 3; % 0: make joint vectors; 1:GMM training; 2:S2H; 3:H2S; 4:PCA; 5:graph; 6:sample for demo

% PCA info
EigenParamDirH = 'J:\!gesture\transitionAmong16of28\EigenParam\1';
EigenParamDirS = 'J:\f0generation\!speech\ATR_A\EigenParam_ATR_A';

% load PCA info
[EVecH, EValH, EuH] = loadEigenParam(EigenParamDirH);
[EVecS, EValS, EuS] = loadEigenParam(EigenParamDirS);
clear EigenParamDirH EigenParamDirS


%% make joint vectors
if mode == 0
    disp('--- make joint vectors ---')
    
    %% definition
    % gesture and speech for convert model
    dirH = 'J:\!gesture\transitionAmong16of28\dgvs';
    %dirS = 'J:\!speech\Japanese5vowels\isolated\suzuki\16k\scep18';
    dirS = 'J:\F0generation\!speech\vowels\feature_scep-f0-vs';
    %dirH = '/Volumes/RESEARCH/_gesture/transitionAmong16of28/dgvs';
    %dirS = '/Volumes/RESEARCH/_speech/Japanese5vowels/isolated/suzuki/16k/scep18';

    consonant = ['b', 'm', 'n', 'p', 'r'];
    vowel = ['a', 'i', 'u', 'e', 'o'];
 
    % the directory for the output files
    %dirOut = '/Volumes/RESEARCH/ProbabilisticIntegrationModel/distortion/bestModel/1113/vowels';
    dirOut = 'J:\f0generation\preliminaryTest\1113_withPCA\joint11_HS';

    % log
    if ismac == 1
    fname_log = [dirOut '/log.txt'];
    else
    fname_log = [dirOut '\log.txt'];
    end

    flog  = fopen(fname_log, 'wt');

    % sample num = 3774
    %a = '28';
    %i = '02';
    %u = '27';
    %e = '22';
    %o = '14';

    % sample num = 1113
    a = '28';
    i = '07';
    u = '01';
    e = '21';
    o = '15';

    fprintf(flog, 'a: %s\t', a);
    fprintf(flog, 'i: %s\t', i);
    fprintf(flog, 'u: %s\t', u);
    fprintf(flog, 'e: %s\t', e);
    fprintf(flog, 'o: %s\t', o);
    fprintf(flog, '\n');
 
    
    %% make joint vectors
    Jnum = 0;
    for nH = 1:1; % the number of dgvs data set
    %for nS = 1:3; % the number of scep data set
    nS = nH;
    % log
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
            %disp([vowel(nV1) vowel(nV2) ': gesture-' num2str(nH) '; speech-' num2str(nS)])
            J = makeJoint(H, S, 1, 0); % with scep0, J = [H, S]
            frameJ = size(J, 2);
            Jnum = Jnum + frameJ;

            % save augmented vector
            if ismac == 1
                fname_joint = [dirOut '/' vowel(nV1) vowel(nV2) '.joint'];
            else
                fname_joint = [dirOut '\' vowel(nV1) vowel(nV2) '.joint'];
            end
            fjoint = fopen(fname_joint, 'wb');
            for ii = 1:frameJ
                fwrite(fjoint, J(:, ii), 'float');
            end
            fclose(fjoint);
    % log
    fprintf(flog, 'joint,%4d\n', frameJ);
            clear fname_joint fjoint frameJ ii;

        end % nV2
    end % nV1
    end % nH and nS

    fprintf(flog, '\ntotal frame number in all joint vectors: %d\n', Jnum);
    fclose(flog);
    clear flog
    clear H S J Jnum
    clear a i u e o vowel
    clear nH nS nV1 nV2


%% GMM training
elseif mode == 1
    H2S = 1; % 0:S2H, 1:H2S
    MIX = 8;
    
    disp('--- GMM training ---')
    Y = loadBinDir('J:\f0generation\preliminaryTest\1113_withPCA\joint11_HS', 'float', 36);
    Y = Y';
    %[EVecH, EValH, EuH] = loadEigenParam(EigenParamDirH);
    %[EVecS, EValS, EuS] = loadEigenParam(EigenParamDirS);

    if H2S == 1
        % for HS joint data
        YH = Y(:, 1:18);
        YS = Y(:, 19:36);
    else
        % for SH joint data
        YH = Y(:, 19:36);
        YS = Y(:, 1:18);
    end
    
    YH = PCA_Trans(YH, EVecH, EuH, 18);
    YS = PCA_Trans(YS, EVecS, EuS, 18);

    if H2S == 1
        % for HS joint data
        Y = [YH, YS];
    else
        % for SH joint data
        Y = [YS, YH];
    end
    clear YH YS
    
    disp(['the number of frames in training data: ' num2str(size(Y, 1))])
    disp(['the number of mixtures: ' num2str(MIX)])
    
    tic
    objH2Smodel = trainGMM(Y, MIX, 0); % full
    save('J:\f0generation\preliminaryTest\1113_withPCA\objH2Smodel-8', 'objH2Smodel');
    toc
    
    clear Y mode MIX
    

%% S2H
elseif mode == 2
    disp('--- S2H conversion ---')
    
    %% definition
    % gesture/speaker model
    dirInModel = 'J:\f0generation';
    %dirInModel = '/Volumes/RESEARCH/ProbabilisticIntegrationModel/distortion';

    % S2Hmodel
    load('J:\f0generation\preliminaryTest\1113_withPCA\objH2Smodel-8');

    %% gesture and speech for convert model
    dirVowels     = 'J:\f0generation\!speech\vowels\feature_scep-f0-vs\1';
    %dirConsonants = 'J:\!speech\JapaneseConsonants\n\suzuki\16k\scep18\1';
    dirConsonants_ = 'J:\f0generation\!speech';

    % the directory for the output files
    dirOutSynDgv = 'J:\f0generation\preliminaryTest\1113_withPCA\synDgv';

    % % input speech
    % dirConsonant_ = 'I:\ProbabilisticIntegrationModel\distortion\scep';
    vowel = ['a', 'i', 'u', 'e', 'o'];
    consonant = ['b', 'm', 'n', 'p', 'r'];

    % Saito's method
    alpha = 1; % weight factor for speaker model
    it = 3; % number or iteration
    updatemethod = 1; %0- using target responsibility 1- using joint responsibility
  
    SAMPLING_FREQ = 1; % assumed sampling frequency of DataGlove
    
    % % PCA parameters
    % [EVec, EVal, Eu] = loadEigenParam(EigenParamDir);


    %% S2H
    % Gesture model
    % according to the preliminary test, the optimal mixture number is 64
    % covariance of objGestureModel is diagonal
    if ismac == 1
        %load([dirInModel '/objGestureModel-64']);
        load([dirInModel '/objGestureModel-32_withPCA']);
    else
        %load([dirInModel '\objGestureModel-64']);
        load([dirInModel '\objGestureModel-32_withPCA']);
    end
    
    nS = 1;
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
                            fname_dgvSyn = [dirOutSynDgv '/' fname '.dgv'];
                            fname_dgvLog = [dirOutSynDgv '/' fname '.txt'];
                        else
                            fname_scepIn = [dirConsonants '\' fname '.feature'];
                            fname_dgvSyn = [dirOutSynDgv '\' fname '.dgv'];
                            fname_dgvLog = [dirOutSynDgv '\' fname '.txt'];
                        end

                    else
                        mora  = sprintf('%s%s', vowel(nV1), vowel(nV2));
                        fname = mora;
                        if ismac == 1
                            fname_scepIn = [dirVowels '/' fname '.feature'];
                            fname_dgvSyn = [dirOutSynDgv '/' fname '.dgv'];
                            fname_dgvLog = [dirOutSynDgv '/' fname '.txt'];
                        else
                            fname_scepIn = [dirVowels '\' fname '.feature'];
                            fname_dgvSyn = [dirOutSynDgv '\' fname '.dgv'];
                            fname_dgvLog = [dirOutSynDgv '\' fname '.txt'];
                        end
                    end
                    disp(mora)
        
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

    
%% H2S    
elseif mode == 3
    disp('--- H2S conversion ---')
    
    %% definition
    % dirH = 'J:\_gesture\transitionAmong16of28\dgvs\3';
    % % dirS = 'J:\_speech\Japanese5vowels\isolated\suzuki\16k\scep18';
    % % %dirH = '/Volumes/RESEARCH/_gesture/transitionAmong16of28/dgvs';
    % % %dirS = '/Volumes/RESEARCH/_speech/Japanese5vowels/isolated/suzuki/16k/scep18';

    % speaker model
    dirInModel = 'J:\f0generation\preliminaryTest\1113_withPCA';
    %dirInModel = '/Volumes/RESEARCH/ProbabilisticIntegrationModel/distortion';

    % Speaker model
    % according to the preliminary test, the optimal mixture number is 128
    % covariance of objSpeakerModel is diagonal
    if ismac == 1
        load([dirInModel '/objSpeakerModel-128']);    
    else
        load([dirInModel '\objSpeakerModel-128']);
    end

    % H2Smodel
    %load('J:\F0generation\withPCA\objS2Hmodel-8');
    load('J:\f0generation\preliminaryTest\1113_withPCA\objH2Smodel-8');    

    % Saito's method
    alpha = 1; % weight factor for speaker model
    it = 3; % number or iteration
    updatemethod = 1; %0- using target responsibility 1- using joint responsibility

    %ENR = 2.5; % the energy of synthesized speech
    
    WINDOW = 40;  % window size for smoothing
    THRES  = 0.0; % threshold for vocing strength

    dirIn  = 'J:\f0generation\preliminaryTest\1113_withPCA\synDgv_demo';
    dirOut = 'J:\f0generation\preliminaryTest\1113_withPCA\synFeature_thres75_demo';

    
    %% H2S
    dirlist = dir([dirIn '\*.dgv']);
    dirlength = length(dirlist);

    ii = 1;
    while ii < dirlength + 1
        % except ".", "..", "DS_Store"
        if length(dirlist(ii).name) > 3 
            filename = strrep(dirlist(ii).name, '.dgv', '');
            disp(filename)

            if ismac == 1
                fname_dgvs      = [dirIn '/' filename '.dgv'];
                fname_scep      = [dirOut '/' filename '.feature'];
                fname_scepLog   = [dirOut '/' filename '.txt'];
                fname_wav       = [dirOut '/' filename '.wav'];
            else
                fname_dgvs      = [dirIn '\' filename '.dgv'];
                fname_scep      = [dirOut '\' filename '.feature'];
                fname_scepLog   = [dirOut '\' filename '.txt'];
                fname_wav       = [dirOut '\' filename '.wav'];
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
            scep_ = scep_';
            
            %scep  = conv2scep(scep_, ENR);
            %scep  = smoothScep(scep, 19, 30);
            scep = scep_;
            clear scep_

            % write out
            frameScep = size(scep, 2);
            fod = fopen(fname_scep, 'wb');
            j = 1;
            while j < frameScep + 1
                fwrite(fod, scep(:, j), 'float');
                j = j + 1;
            end
            clear scep
            fclose(fod);
            
            % feature to wav
            feature2wav(fname_scep, WINDOW, THRES);
        end
        ii = ii + 1;
    end % while
    clear it alpha input j
    clear dirIn dirInModel dirOut fod
    clear filename fname_dgvs fname_scep fname_scepLog fname_wav frameScep 

    
%% PCA
elseif mode == 4
    X = loadBinDir('J:\F0generation\!speech\ATR_A\feature_scep-f0-vs','float',18);
    X = X';
    getEigenParam(X, 'J:\F0generation\withPCA\EigenParam_ATR_A');

    
%% draw graph
elseif mode == 5
    clear EValH EVecH EuH 
    clear EValS EVecS EuS
    
    % vowels
    %dirTgt = 'J:\F0generation\!speech\vowels';
    % consonant
    %dirTgt = 'J:\F0generation\!speech\p';
    % ATR
    dirTgt = 'J:\f0generation\!speech\ATR_J';

    %dirEst = 'J:\F0generation\bestModel\1113\synFeature';
    dirEst = 'J:\f0generation\preliminaryTest\1113_withPCA\synFeature_thres75_demo';
    mora = 'j01';

    %% make filenames
    % target
    % vowel/consonant
    %fF0_tgt = [dirTgt '\f0\1\' mora '.f0'];
    %fVS_tgt = [dirTgt '\vs_usingNCCF\1\' mora '.vs'];
    % ATR
    fF0_tgt = [dirTgt '\f0\' mora '.f0'];
    fVS_tgt = [dirTgt '\vs_usingNCCF\' mora '.vs'];

    clear dirTgt

    % estimated
    fF0_est = [dirEst '\' mora '.f0'];
    fVS_est = [dirEst '\' mora '.vs'];
    clear dirEst


    %% load files

    % target
    F0_tgt = load(fF0_tgt);
    VS_tgt = load(fVS_tgt);
    clear fF0_tgt fVS_tgt

    % estimated
    F0_est = load(fF0_est);
    VS_est = load(fVS_est);
    clear fF0_est fVS_est


    %% arrange data format
    CUT  = 40;
    fmaxF0 = length(F0_tgt);
    fmaxVS = length(VS_tgt);
    if fmaxF0 ~= fmaxVS
        F0_tgt(fmaxF0) = [];
    end
    fmax = fmaxVS;
    clear fmaxF0 fmaxVS

    F0_tgt(fmax-CUT+1:fmax) = [];
    F0_tgt(1:CUT) = [];
    VS_tgt(fmax-CUT+1:fmax) = [];
    VS_tgt(1:CUT) = [];
    clear CUT fmax

    hold on
        %plot(F0_tgt, 'r')
        %plot(F0_est, 'b')
        plot(VS_tgt, 'r')
        plot(VS_est, 'b')
    hold off

    
%% sample for demo
elseif mode == 6
    %% definition
    dirInS = 'J:\f0generation\!speech\ATR_J\feature_scep-f0-vs';
    dirSynDgv = 'J:\f0generation\preliminaryTest\1113_withPCA\synDgv_demo';
    fname = 'j03';
    
    % gesture/speaker model
    dirInModel = 'J:\f0generation';
    %dirInModel = '/Volumes/RESEARCH/ProbabilisticIntegrationModel/distortion';

    % S2Hmodel
    load('J:\f0generation\preliminaryTest\1113_withPCA\objH2Smodel-8');
    % Saito's method
    alpha = 1; % weight factor for speaker model
    it = 3; % number or iteration
    updatemethod = 1; %0- using target responsibility 1- using joint responsibility
  
    SAMPLING_FREQ = 1; % assumed sampling frequency of DataGlove
    
    
    %% S2H
    % Gesture model
    % according to the preliminary test, the optimal mixture number is 64
    % covariance of objGestureModel is diagonal
    if ismac == 1
        %load([dirInModel '/objGestureModel-64']);
        load([dirInModel '/objGestureModel-32_withPCA']);
    else
        %load([dirInModel '\objGestureModel-64']);
        load([dirInModel '\objGestureModel-32_withPCA']);
    end
    
    fname_scepIn = [dirInS '\' fname '.feature'];
    fname_dgvSyn = [dirSynDgv '\' fname '.dgv'];
    fname_dgvLog = [dirSynDgv '\' fname '.txt'];

    % conversion
    input = loadBin(fname_scepIn, 'float', 18);
    input = input';
    input = PCA_Trans(input, EVecS, EuS);
    input = input';

    dgv_ = spkmodel_vc2_(input, objH2Smodel, objGestureModelWithPCA, alpha, it, updatemethod, fname_dgvLog);
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
end