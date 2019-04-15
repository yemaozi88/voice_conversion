%
% 2010-04-22
% check recorded dgv files
%

clear all, clc, fclose('all');

% nn = 'n27';
% cc = 'no';
% dirIn = 'C:';
% dirOut_ = 'C:\hmts\100430_INTERSPEECH';
% 
% 
% %% check mean
% 
% finA1 = [dirIn '\1\' cc '.dgv'];
% finA2 = [dirIn '\2\' cc '.dgv'];
% finA3 = [dirIn '\3\' cc '.dgv'];
% 
% A1 = loadBin(finA1, 'uchar', 26);
% A2 = loadBin(finA2, 'uchar', 26);
% A3 = loadBin(finA3, 'uchar', 26);
% 
% %B1 = loadBin(['C:\hmts\100419_ProgressReport\' nn '\dgv\n\1\' cc '.dgv'], 'uchar', 26);
% %B2 = loadBin(['C:\hmts\100419_ProgressReport\' nn '\dgv\n\2\' cc '.dgv'], 'uchar', 26);
% %B3 = loadBin(['C:\hmts\100419_ProgressReport\' nn '\dgv\n\3\' cc '.dgv'], 'uchar', 26);
% 
% [mean(A1')', mean(A2')', mean(A3')']
% 
% 
% %% cut the first bad data
% dirOut  = [dirOut_ '\' nn '\dgv\n'];
% 
% foutA1 = [dirOut '\1\' cc '.dgv'];
% foutA2 = [dirOut '\2\' cc '.dgv'];
% foutA3 = [dirOut '\3\' cc '.dgv'];
% 
% cut1dgv(finA1, foutA1);
% cut1dgv(finA2, foutA2);
% cut1dgv(finA3, foutA3);
% 
% 
% %% display STM
% STM1 = getSTM(finA1);
% STM2 = getSTM(finA2);
% STM3 = getSTM(finA3);
% 
% fmax1 = size(STM1, 2);
% fmax2 = size(STM2, 2);
% fmax3 = size(STM3, 2);
% 
% fmax = max([fmax1, fmax2, fmax3]);
% STM1_ = zeros(1, fmax);
% for ii = 1:fmax1 
%     STM1_(1, ii) = STM1(1, ii);
% end
% STM2_ = zeros(1, fmax);
% for ii = 1:fmax2 
%     STM2_(1, ii) = STM2(1, ii);
% end
% STM3_ = zeros(1, fmax);
% for ii = 1:fmax3 
%     STM3_(1, ii) = STM3(1, ii);
% end
% 
% hold on
% plot(STM1_, 'r')
% plot(STM2_, 'g')
% plot(STM3_, 'b')
% hold off
% 
% STM = STM1_ + STM2_ + STM3_;
% %plot(STM)

% au(28-25), ui(25-07)
%y=loadBin([dirOutJoint '/au_11.joint'], 'float', 36);
%mean(y')'

% dirH = '/Users/kunikoshi/research/gesture/transitionAmong16of28/dgvs';
% dirS = '/Users/kunikoshi/research/speech/Japanese5vowels/isolated/suzuki/16k/scep18';
% datasetH = '1';
% datasetS = '1';
% 
% datanumH = '28-25';
% datanumS = 'au';
% 
% filenameH = [dirH '/' datasetH '/' datanumH '.dgvs'];
% filenameS = [dirS '/' datasetS '/' datanumS '.scep'];
% 
% H = loadBin(filenameH, 'uchar', 26);
% %H = H(5:22, :);
% S = loadBin(filenameS, 'float', 19);
% %S = S(2:19, :);
% 
% J = makeJoint(H, S, 1);


%% graph for cepstrum distortion
% mora = 'aa.scep';
% dirS = 'C:\research\speech\Japanese5vowels\isolated\suzuki\16k\scep18\1';
% dirInOrg = [dirS '\' mora];
% 
% dirSyn = 'C:\research\ProbabilisticIntegrationModel\S2H-H2S_ERRV20_ERRC20_thres5_mix32\28-07-04-13-14\synScep';
% dirInSyn = [dirSyn '\' mora];
% 
% ORG = loadBin(dirInOrg, 'float', 19);
% SYN = loadBin(dirInSyn, 'float', 19);
% 
% deg = 9;
% hold on
%     plot(ORG(deg, 1:100)', 'b');
%     plot(SYN(deg, 1:100)', 'r');
% hold off


%% the comparison between GMM for H2S and S2H
%dirIn = 'C:\research\ProbabilisticIntegrationModel\S2H-H2S_distortion\preliminaryTest\28-07-04-13-14';
%dirIn = 'C:\research\ProbabilisticIntegrationModel\S2H-H2S_distortion\5840_28-22-11-07-21';

% MIX = 8;
% 
% H2S = loadBinDir([dirIn '\joint_H2S\H1S1'], 'float', 36);
% objH2Smodel = trainGMM(H2S', MIX, 0);
% save([dirIn '\objH2Smodel_mix' num2str(MIX)], 'objH2Smodel');
%  
% S2H = loadBinDir([dirIn '\joint_S2H\H1S1'], 'float', 36);
% objS2Hmodel = trainGMM(S2H', MIX, 0);
% save([dirIn '\objS2Hmodel_mix' num2str(MIX)], 'objS2Hmodel');

% dirS = 'C:\research\speech\Japanese5vowels\isolated\suzuki\16k\scep18';
% dirConsonant_ = 'C:\research\speech\JapaneseConsonants';
% consonant = ['b', 'm', 'n', 'p', 'r'];
% vowel = ['a', 'i', 'u', 'e', 'o'];
% 
% % Saito's method
% alpha = 1; % weight factor for speaker model
% it = 3; % number or iteration
% updatemethod = 1; %0- using target responsibility 1- using joint responsibility
% 
% ENR = 2.5; % the energy of synthesized speech
% SAMPLING_FREQ = 1; % assumed sampling frequency of DataGlove
% 
% load('C:\research\ProbabilisticIntegrationModel\S2H-H2S_distortion\objGestureModel-64');
% load('C:\research\ProbabilisticIntegrationModel\S2H-H2S_distortion\objSpeakerModel-64');
% %load('C:\research\ProbabilisticIntegrationModel\S2H-H2S_distortion\preliminaryTest\28-07-04-13-14\objS2Hmodel_mix8');
% load('C:\research\ProbabilisticIntegrationModel\S2H-H2S_distortion\5840_28-22-11-07-21\objS2Hmodel');
% 
% dirOutSynDgv  = 'C:\research\ProbabilisticIntegrationModel\S2H-H2S_distortion\5840_28-22-11-07-21\synDgv';
% dirOutSynScep = 'C:\research\ProbabilisticIntegrationModel\S2H-H2S_distortion\5840_28-22-11-07-21\synScep';
% 
% % consonant
% dirConsonant_ = 'C:\research\speech\JapaneseConsonants';
% consonant = ['b', 'm', 'n', 'p', 'r'];
% 
% %nH = 1;
% nS = 1;
% for nC = 1:5; % n
% for nV = 1:5
% % S2H
%     %mora = sprintf('%s%s', vowel(nV), vowel(nV));
%     mora = sprintf('%s%s', consonant(nC), vowel(nV));
% disp(['---' mora ' ---'])
% 
%     if ismac == 1
%         fname_scepIn = [dirConsonant_ '/' consonant(nC) '/suzuki/16k/scep18/' num2str(nS) '/' mora '.scep'];
%     else
%         fname_scepIn = [dirConsonant_ '\' consonant(nC) '\suzuki\16k\scep18\' num2str(nS) '\' mora '.scep'];
%     end
%                 
%     if ismac == 1
%         %fname_scepIn    = [dirS '/' num2str(nS) '/' mora '.scep'];
%         fname_dgvSyn    = [dirOutSynDgv '/' mora num2str(nS) '.dgv'];
%         fname_dgvLog    = [dirOutSynDgv '/' mora num2str(nS) '.txt'];
%         fname_scepSyn   = [dirOutSynScep '/' mora num2str(nS) '.scep'];
%         fname_scepLog   = [dirOutSynScep '/' mora num2str(nS) '.txt'];
%         fname_wav       = [dirOutSynScep '/' mora num2str(nS) '.wav'];
%     else
%         %fname_scepIn    = [dirS '\' num2str(nS) '\' mora '.scep'];
%         fname_dgvSyn    = [dirOutSynDgv '\' mora num2str(nS) '.dgv'];
%         fname_dgvLog    = [dirOutSynDgv '\' mora num2str(nS) '.txt'];
%         fname_scepSyn   = [dirOutSynScep '\' mora num2str(nS) '.scep'];
%         fname_scepLog   = [dirOutSynScep '\' mora num2str(nS) '.txt'];
%         fname_wav       = [dirOutSynScep '\' mora num2str(nS) '.wav'];
%     end
%     
%     input = loadBin(fname_scepIn, 'float', 19);
%     input = input(2:19, :); % remove energy
%     dgv_ = spkmodel_vc2(input, objS2Hmodel, objGestureModel, alpha, it, updatemethod, fname_dgvLog);
%     % dgv_ holds all results at every step
%     dgv_ = dgv_{1, it};
%     dgv  = conv2dgv(dgv_, SAMPLING_FREQ);
%     frameDgv = size(dgv, 2);
%     
%     fout = fopen(fname_dgvSyn, 'wb');
%     for ii = 1:frameDgv
%         fwrite(fout, dgv(:, ii), 'uchar');
%     end
%     fclose(fout);
% disp('S2H completed');
%     
% % H2S
%     dgv = dgv(5:22, :); % remove energy
%     
%     scep_ = spkmodel_vc2_(dgv, objS2Hmodel, objSpeakerModel, alpha, it, updatemethod, fname_scepLog);
%     % scep_ holds all results at every step
%     scep_ = scep_{1, it};
%     scep  = conv2scep(scep_, ENR);
%     scep = smoothScep(scep, 19, 30);
%     clear scep_
%    
%     % write out scep
%     frameScep = size(scep, 2);
%     fod = fopen(fname_scepSyn, 'wb');
%     j = 1;
%     while j < frameScep + 1
%         fwrite(fod, scep(:, j), 'float');
%         j = j + 1;
%     end
%     fclose(fod);
% 
%     % write out wav synthesized using scep
%     scep2wav(scep, fname_wav, 'suzuki'); 
% disp('H2S completed');
% 
% % cepstral distortion
%     NMmean = distortion2(fname_scepIn, fname_scepSyn);
% disp(['cepstral distortion: ' num2str(NMmean)]);
% end %nV
% end %nC

