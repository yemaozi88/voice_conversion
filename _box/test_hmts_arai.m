% 2008-8-22
% test realtime hmts project

clear all

%- deg of jointvector
NUM = 22; % number of signal from dataglove
SNS = 18; % number of sensor
DEG = 17; % degree of cepstrum, usually SNS - 1
JNT = SNS + DEG + 1; % degree per a joint vector


%% load matrix data
cd 'H:\u-tokyo\HMTS-win\kuni_model_ASJ\speech\suzuki\16k\joint\matrix1'
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


%% get cepstrum vector from input
%input = importdata('C:\Users\kunikoshi\Documents\HMTS-win\demo\dgv\test\ao.txt');
%fnum = length(input(:,1));

finM = fopen('H:\u-tokyo\HMTS-win\kuni_model_ASJ\hand\ai_open\ai2.dgv','rb');

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
            B(:, s) = B_(5:SNS+4, step);
        end
    elseif t == fmax;
        for s = step_buf:step
            B(:, s) = (B_(5:SNS+4, t) - B_(5:SNS+4, t-1)) * (s-step_buf)/B_(1, t) + B_(5:SNS+4, t-1);
        end
    else
        for s = step_buf:step-1
            B(:, s) = (B_(5:SNS+4, t) - B_(5:SNS+4, t-1)) * (s-step_buf)/B_(1, t) + B_(5:SNS+4, t-1);
        end
    end
    step_buf = step;
end

cep = zeros(DEG+1, total_time);
for t = 1:total_time;
    f = mean_w + covar_mul * (B(:, t) - mean_v);
    cep(:, t) = f;
end


%% resynthesis speech
% get spectrogram from cepstrum time series
sgram = cep2sgram(cep);


% define sampling rate
fs = 16000; % 10[ms] = 100[Hz]


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
wavwrite(sy/32768, fs, 16, 'H:\u-tokyo\HMTS-win\kuni_model_ASJ\speech\suzuki\16k\syn\ai2.wav');