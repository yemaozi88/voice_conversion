% 2008-4-28
% make gaussian model
% 2008-7-11
% apply for dgv data

clear all

%- deg of jointvector
NUM = 22; % number of signal from dataglove
SNS = 18; % number of sensor
DEG = 17; % degree of cepstrum, usually SNS - 1
JNT = 36; % degree per a joint vector


%% load joint vector
% %cd 'C:\Documents and Settings\kunikoshi\My Documents\HMTS-win\080712kuni_model2';
% cd 'C:\Users\kunikoshi\Documents\voice_data\suzuki'
% dirlist = dir('joint_16k');
% dirlength = length(dirlist);
% 
% A = [];
% 
% i = 1;
% while i < dirlength + 1
% 	%% except ".", "..", "DS_Store"
%   	if length(dirlist(i).name) > 3 
%         filename = dirlist(i).name;
%    		finJ = fopen(['joint_16k\' filename], 'rb');
%         A_ = [fread(finJ, [JNT, inf], 'float')];
%         A = [A, A_];
%     end
% 	i = i + 1;
% end


%% get mean vector and covariance matrix
% m = mean(A')';
% C = cov(A');
% iC = inv(C);


%% load matrix data
cd 'C:\Users\kunikoshi\Documents\HMTS-win\kuni_model_2\speech\suzuki\16k\joint\matrix2'
mean_v = importdata('mean_v.txt');
mean_w = importdata('mean_w.txt');
covar_mul_buf = importdata('covar_mul.txt');

% set covar_mul_buf value to covar_mul 
covar_mul = zeros(18, 18);
for(r = 1:SNS)
    for(c = 1:SNS)
        covar_mul(c, r) = covar_mul_buf((r - 1)*SNS + c, 1); 
    end
end


%% load input dgv file
%finM = fopen('C:\Documents and Settings\kunikoshi\My Documents\HMTS-win\080712kuni_model2\dgv\ai.dgv','rb');
%finM = fopen('C:\Users\kunikoshi\Documents\HMTS-win\080712kuni_model2\ai_test\ai5.dgv','rb');
finM = fopen('C:\Users\kunikoshi\Documents\HMTS-win\kuni_model_2\hand\ai_open\ai2.dgv','rb');

B_ = fread(finM, [NUM + 4, inf], 'uchar');
fmax = length(B_(1,:));
T = sum(B_');
total_time = T(1);
 
B = zeros(SNS, total_time);
 
step = 0;
step_buf = 0;
for t = 1:fmax;
    step = step_buf + B_(1, t);
    if t == 1; % step_buf == 0;
        for s = 1:step-1
            B(:, s) = B_(5:22, step);
        end
    elseif t == fmax;
        for s = step_buf:step
            B(:, s) = (B_(5:22, t) - B_(5:22, t-1)) * (s-step_buf)/B_(1, t) + B_(5:22, t-1);
        end
    else
        for s = step_buf:step-1
            B(:, s) = (B_(5:22, t) - B_(5:22, t-1)) * (s-step_buf)/B_(1, t) + B_(5:22, t-1);
        end
    end
    step_buf = step;
end


%% get cepstrum vector from B
scep = zeros(DEG+1, total_time);
for t = 1:total_time;
    %f = m(SNS+1:JNT) + C(SNS+1:JNT, 1:SNS)/C(1:SNS, 1:SNS) * (B(:, t) - m(1:SNS));
    f = mean_w + covar_mul * (B(:, t) - mean_v);
    scep(:, t) = f;
end

 
%% resynthesis speech
% get spectrogram from cepstrum time series
sgram = cep2sgram(scep);
 
% % original...
% %[x, fs] = wavread('C:\Users\kunikoshi\Documents\matsuura\ai.wav');
% %[f0raw, ap, prmF0] = exstraightsource(x, fs);
 
% define sampling rate
fs = 16000;

% make ap matrix
ap2 = sgram;
% for t = 1:total_time
%     if t > length(ap)
%         ap2(:, t) = ap(:, length(ap));
%     else
%         ap2(:, t) = ap(:, t);
%     end
% %     for s = 1:513
% %         ap(s, t) = 0;
% %     end
% end
for t = 1:total_time
    for s = 1:513
        ap2(s, t) = -20;
    end
end

% make f0 vector
f0 = zeros(1, total_time);
for t = 1:total_time
     f0(1, t) = 140;
end


%% write to wav file
sy = exstraightsynth(f0, sgram, ap2, fs);
wavwrite(sy/32768, fs, 16, 'C:\Users\kunikoshi\Documents\HMTS-win\kuni_model_2\speech\suzuki\48k\syn\16k_ai2.wav');

%fclose(finJ);
fclose(finM);