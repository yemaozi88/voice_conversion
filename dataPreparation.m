%
% 2009/05/02
% do processes to all files in the directory
% 
% AUTHOR
% Aki Kunikoshi (D1)
% yemaozi88@gmail.com
%

fclose all, clear all, clc

%% definition
default;
%dirIn  = 'd:\akunikoshi\Research\MRI\speechWithoutNoise\cmu\';
%dirOut = 'd:\akunikoshi\Research\MRI\speechWithoutNoise\cmu_resyn\';

% extract features
%dirIn_       = 'J:\F0generation\withPCA\synFeature';
%dirScep    = 'J:\ProbabilisticIntegrationModel\F0generation\ATR_J\scep18';
%dirResyn  = 'J:\!speech\ATR503\suzuki\16k_from_48k\resyn18';
%dirVS       = 'J:\ProbabilisticIntegrationModel\F0generation\ATR_J\vs_usingNCCF';
%dirF0       = 'J:\ProbabilisticIntegrationModel\F0generation\ATR_A\f0_usingNCCF';
%dirIF0     = 'J:\ProbabilisticIntegrationModel\F0generation\ATR_J\f0interp';
%dirOut     = 'J:\ProbabilisticIntegrationModel\F0generation\ATR_J\feature_scep-f0-vs';

% DTW
% dirInX   = 'J:\VoiceConversion\STRAIGHT-based\fws\scep18+delta_train_reduced1of5';
% dirInY   = 'J:\VoiceConversion\STRAIGHT-based\mht\scep18+delta_train_reduced1of5';
% dirCSV   = 'J:\VoiceConversion\STRAIGHT-based\fws2mht_withDelta_reduced1of5\csv_train';
% dirJoint = 'J:\VoiceConversion\STRAIGHT-based\fws2mht_withDelta_reduced1of5\joint_train';

% for cc = ['m', 'n'];
%     tic
%     for nn = 6:10;
% 
% dirIn = [dirIn_ '\' cc '\suzuki\16k'];
%dirIn = dirIn_;
%    nnStr = num2str(nn);
%dirWav  = [dirIn '\' cc '-cmv\wav\' num2str(nn)];
%dirScep = [dirIn '\' cc '-cmv\scep18\' num2str(nn)];
% dirWav   = [dirIn '\wav\' num2str(nn)];
% dirF0    = [dirIn '\f0_\' num2str(nn)];
% dirScep  = [dirIn '\scep16\' num2str(nn)];
% dirResyn = [dirIn '\resyn16\' num2str(nn)];
% dirWav   = [dirIn '\wav\' nnStr];
% dirF0    = [dirIn '\f0\' nnStr];
% dirScep  = [dirIn '\scep16\' nnStr];
% dirResyn = [dirIn '\resyn16\' nnStr];

%dirVS   = [dirIn '\vs_usingNCCF\'];
%dirF0   = [dirIn '\f0_usingNCCF\'];
%dirIF0  = [dirIn '\f0interp\' num2str(nn)];
%dirOut  = [dirIn '\resyn16\' num2str(nn)];


%% HMTS

%dirModel = ['C:\hmts\100430_INTERSPEECH\' NN '\matrix_' num2str(MIX)'];
%dirModel = 'C:\hmts\100430_INTERSPEECH\vowels\matrix_2_PCA-8_reduced-10';
%Model = loadGMM_Qiao(dirModel, 8, 18, 2, 5);
%[Evec, Eval, u] = loadEval(['C:\hmts\100430_INTERSPEECH\' NN '\Eval']);
%[Evec, Eval, u] = loadEval('C:\hmts\100430_INTERSPEECH\n16\Eval');
%EigenParamDirH = 'J:\!gesture\transitionAmong16of28\EigenParam\all';
%[EVecH, EValH, uH] = loadEigenParam(EigenParamDirH);
%clear EigenParamDirH EValH

% load('J:\ProbabilisticIntegrationModel\distortion\objGestureModel-64');
% load('J:\ProbabilisticIntegrationModel\distortion\5840\objS2Hmodel');
% load('J:\ProbabilisticIntegrationModel\distortion\5840\joint\obj\jointModel_mix64_obj');
% 
% % Saito's method
% alpha = 1; % weight factor for speaker model
% it = 3; % number or iteration
% updatemethod = 1; %0- using target responsibility 1- using joint responsibility
% 
% ENR = 2.5; % the energy of synthesized speech
% SAMPLING_FREQ = 1; % assumed sampling frequency of DataGlove


%% f0 generation
%dirIn  = 'J:\F0generation\withPCA\synFeature';
%WINDOW = 40;  % window size for smoothing
%THRES  = 0.75; % threshold for vocing strength


%% directory processing
dirlist = dir([dirIn_wav '\*.wav']);
for ii = 1:length(dirlist)
%for ii = 10:10
%for ii = 1:3 % test data is the first one
	% except ".", "..", "DS_Store"
  	if length(dirlist(ii).name) > 3
        %filename = dirlist(ii).name;
        fin_name  = strrep(dirlist(ii).name, '.wav', '');
        fout_name = strrep(fin_name, 'fws', 'mht');
        %[pathstr, name, ext] = fileparts(ffeature);
        %disp(['No. ' num2str(ii) ' : ' filename])
        %disp([num2str(nn) ' - ' filename])
        
        fin_wav  = [dirIn '\' fin_name '.wav'];
        fout_wav = [dirOut '\' fout_name '.wav'];
        %disp([num2str(nn) ':' filename])
        fprintf('%s --> %s\n', fin_name, fout_name);
        %disp(fout)
        
        %% extract features
      %fin  = [dirWav '\' filename '.wav'];
      %fout = [dirScep '\' filename '.scep'];
%         filename = strrep(dirlist(ii).name, 'fws_m_', '');
%         filename = strrep(filename, '_0', '');
%         filename = strrep(filename, '.wav', '');
        %fWav   = [dirWav '\' filename '.wav'];
        %fF0    = [dirF0 '\' filename '.f0'];
        %fIF0   = [dirIF0 '\' filename '.if0'];
        %fScep  = [dirScep '\' filename '.scep'];
        %fResyn = [dirResyn '\' filename '.wav'];
        %fVS    = [dirVS '\' filename '.vs'];
        %fout   = [dirOut '\' filename '.wav'];
        
%         disp(fScep)
%         disp(fF0)
%         disp(fVS)
%         disp(fWav)
%         disp(fResyn)

      %% DTW
%         finX        = [dirInX '\' filename '.feature'];
%         finY        = [dirInY '\' filename '.feature'];
%         foutCSV     = [dirCSV '\' filename '.csv'];
%         foutJoint   = [dirJoint '\' filename '.joint'];        
%  
%         X = loadBin(finX, 'float', 36);
%         Y = loadBin(finY, 'float', 36);
        
        % remove energy
        %X = X(2:19, :);
        %Y = Y(2:19, :);
        
        
        %% HMTS
%         input = loadBin(fin, 'float', 19);
%         input = input(2:19, :); % remove energy
%         dgv_ = spkmodel_vc2(input, objS2Hmodel, objGestureModel, alpha, it, updatemethod, fname_dgvLog);
%         % dgv_ holds all results at every step
%         dgv_ = dgv_{1, it}; % dgv_ is PCAed dgv
%         %dgvPCA = dgv_
%         
%         dgv  = conv2dgv(dgv_, SAMPLING_FREQ);
%         frameDgv = size(dgv, 2);
% 
%         fout_ = fopen(fout, 'wb');
%         for ii = 1:frameDgv
%             fwrite(fout_, dgv(:, ii), 'uchar');
%         end
%         fclose(fout_);
%         clear filename fname
%         clear input dgv_ dgv frameDgv

        %disp(['fin:' fin])
        %disp(['fout:' fout])

        %% do
        %extractScep(fin, fout, 20);
        %extractFeature(fWav, fF0, fScep, fResyn, 16);
        %extractVSusingNCCF(fWav, fVS, fF0);
        %extractF0(fin, fF0);
        %interpolateF0(fF0, fIF0, 0); % log-scale
        %resyn(fin, fout, 20, 'natural');
        %dtwCal(X, Y, foutCSV);
        %dtwProc(X, Y, 36, foutCSV, foutJoint);
        %mkdata(fScep, 18, fIF0, fVS, fout, 40);
        %mkdata2(fin, 18, fout);
        %addDelta2(fin, 18, 'dgv', fout);
        %HMTS_suzuki(dirMatrix, fin, fout);
        %joint2htk(fin, fout);
        %mappingGMM_Qiao(fin, fout, Model);
        %mappingGMMdelta_Qiao(fin, fout, Model, 5);
        %cut1joint(fin, fout);
        %PCA4joint(fin, fout, Evec, u, 8);
        %reduceData(fin, fout, 'float', 32, 8);
        %mappingPCA2GMM_Qiao(fin, fout, Model, Evec, u, 8);
        %addDelta2joint(fin, fout, 5);
        %feature2wav(fin, WINDOW, THRES);
        %MCDT = distortion(fin, fout);
        %fprintf(fout1, '%s,', filename);
        %fprintf(fout1, '%f\n', MCDT);
        %removeFirstData(fin, fout, 'uchar', 26);
        %extractFeatureWithoutSil(fin, foutF0, foutScep, 18, -4.3, 0);
        %extractStableFrames(fin, fout, EVecH, uH);
    end
end % ii

%     end % nn
%     toc
% end % cc