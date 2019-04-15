%
% 2009-07-10 Morphing
%
% LINK
% STRAIGHT
%
% REMARKS
% Hand Motions are interpolated by every 1ms after loading
% f0 is 140[Hz] constant
% ap is -20 constant
% NUM = 22; % number of signal from dataglove
% SNS = 18; % number of sensor
% DEG = 17; % degree of cepstrum, usually SNS - 1
% JNT = 36; % degree per a joint vector
%

%% morphing object
moA = createMobject;
moB = createMobject;

%% load wav files
[x, fs] = wavread('C:\Documents and Settings\kunikoshi\デスクトップ\rec\aa_44_s.wav');
moA = updateFieldOfMobject(moA, 'waveform', x);
%moA.samplingFrequency = 16000;
[x, fs] = wavread('C:\Documents and Settings\kunikoshi\デスクトップ\rec\ii_44_s.wav');
moB = updateFieldOfMobject(moB, 'waveform', x);
%moB.samplingFrequency = 16000;

%% analyse
moA = executeSTRAIGHTanalysisM(moA);
moB = executeSTRAIGHTanalysisM(moB);

save moA moA
save moB moB

%mObjectdmy = timeFrequencySTRAIGHTmorphing(moA, moB, 0.5, 'log');
%displayMobject(mObjectdmy, 'anchorTimeLocation' ,'Log');
%axis([0 283 0 7000])

%mObjectdmy = timeAlignedDirectSTRAIGHTmorphing(moA, moB, 0.5, 'log');
%displayMobject(mObjectdmy, 'anchorTimeLocation', 'alignedLinSum');
%axis([0 283 0 7000])

%mixedWave = waveformMorphing(moA, moB, 0.5);
%displayMobject(mixedWave,'waveform','mixedWave');

%syneu = executeSTRAIGHTsynthesisM(mixedWave);



%fout = 'mo.wav';
%wavwrite(syneu/32768, syneu, 16, fout);

% %% definition
% NUM = 22; % number of signal from dataglove
% SNS = 18; % number of sensor
% DEG = 17; % degree of cepstrum, usually SNS - 1
% JNT = 36; % degree per a joint vector
% 
% % dr: the directory mean_v, mean_w and covar_mul located
% % fin: a input dgv file from 18 sensor CyberGlove
% % fout: where wav file will be located
% dr = 'J:\hmts\090906_INTERSPEECH\model05\matrix3';
% fin = 'J:\hmts\090906_INTERSPEECH\model05\dgv\open\aiueo_ambiguous1.dgv';
% %fout = 'J:\hmts\090906_INTERSPEECH\model05\syn\aiueo_clear1.wav';
% fout_cep = fopen('J:\hmts\090906_INTERSPEECH\model05\sgram.txt', 'wt');
% 
% 
% %% load matrix data
% fn_mean_v = [dr '\mean_v.txt'];
% fn_mean_w = [dr '\mean_w.txt'];
% fn_covar_mul = [dr '\covar_mul.txt'];
% mean_v = importdata(fn_mean_v);
% mean_w = importdata(fn_mean_w);
% covar_mul_buf = importdata(fn_covar_mul);
% 
% % set covar_mul_buf value to covar_mul 
% covar_mul = zeros(18, 18);
% for(r = 1:SNS)
%     for(c = 1:SNS)
%         covar_mul(c, r) = covar_mul_buf((r - 1)*SNS + c, 1); 
%     end
% end
% 
% 
% %% get cepstrum vector from input  
% finM = fopen(fin,'rb');
%  
% B_ = fread(finM, [NUM + 4, inf], 'uchar');
% fmax = length(B_(1,:));
% T = sum(B_');
% total_time = T(1);
%  
% B = zeros(SNS, total_time);
%  
% step = 0;
% step_buf = 0;
% for t = 1:fmax;
%     step = step_buf + B_(1, t);
%     if t == 1; % step_buf == 0;
%         for s = 1:step-1
%             B(:, s) = B_(5:22, step);
%         end
%     elseif t == fmax;
%         for s = step_buf:step
%             B(:, s) = (B_(5:22, t) - B_(5:22, t-1)) * (s-step_buf)/B_(1, t) + B_(5:22, t-1);
%         end
%     else
%         for s = step_buf:step-1
%             B(:, s) = (B_(5:22, t) - B_(5:22, t-1)) * (s-step_buf)/B_(1, t) + B_(5:22, t-1);
%         end
%     end
%     step_buf = step;
% end
% 
% fclose(finM);
% 
% 
% %% get cepstrum time series
% cep = zeros(DEG+1, total_time);
% for t = 1:total_time;
%     f = mean_w + covar_mul * (B(:, t) - mean_v);
%     cep(:, t) = f;
% end
% 
% 
% %% resynthesis speech
% % get spectrogram from cepstrum time series
% sgram = cep2sgram(cep);
% 
% 
% %% write out cep data into .txt file
% fprintf(fout_cep, '%f\n', sgram(:, 3000));
% 
% % define sampling rate
% %fs = 16000; % 10[ms] = 100[Hz]
% 
% % make ap matrix
% % ap2 = sgram;
% % 
% % for t = 1:total_time
% %     for s = 1:513
% %         ap2(s, t) = -20;
% %    end
% %end
% 
% % make f0 vector
% %f0 = zeros(1, total_time);
% %for t = 1:total_time
% %     f0(1, t) = 140;
% %end
% 
% 
% %% write to wav file
% %sy = exstraightsynth(f0, sgram, ap2, fs);
% %wavwrite(sy/32768, fs, 16, fout);
% 
% 
% %% file close
% fclose(fout_cep);