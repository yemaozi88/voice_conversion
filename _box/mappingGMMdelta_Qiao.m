function mappingGMMdelta_Qiao(fin, fout, Model, L)
% space mapping based on GMM based on Qiao's program
%
% mappingGMM_Qiao(fin, fout, Model)
%
% INPUT
% fin: a input dgv file from 18 sensor CyberGlove
% fout: where wav file will be located
% model: GMM data
% L: delta window size
%
% LINK
% STRAIGHT, cep2sgram
%
% REMARKS
% model can be gotten by getGMM_Qiao or loaded by loadGMM_Qiao.m
%
% Hand Motions are interpolated by every 1ms after loading
% f0 = 140[Hz] constant
% ap = -20 constant
% NUM = 22; % number of signal from dataglove
% SNS = 18; % number of sensor
% DEG = 17; % degree of cepstrum, usually SNS - 1
% JNT = 36; % degree per a joint vector
%
% HISTORY
% 2009-10-25 calculate f(a) part was re-written based on GMM_Mapping.m
% 2010-01-15 functionized
%

%fin = 'C:\hmts\na4.dgv';
%fout = 'C:\hmts\na4.wav';
%L = 5;


%% definition
% this program should be changed using this part
NUM = 22;
SNS = 18;
DEG = 17;
JNT = 54;

%Model = loadGMM_Qiao('C:\hmts\091008_SP\model4-n\matrix_5-32', 18,18, 32,1);


%% get dgv vector from input 
finM = fopen(fin, 'rb');
%finM = fopen('C:\hmts\091008_SP\model4-n\dis_5-32\dgv\n-cmv\3\ni-m.dgv', 'rb');

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
            B(:, s) = B_(5:22, 1);
            %B(:, s) = B_(5:22, step); // it is wrong
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

fclose(finM);

DELTA = dgv2delta(B, 5);
B = [B; DELTA];

% cut both ends
fmax = size(B ,2);
B = B(:, (L+1)+1:fmax-(L+1));


%% calculate the output cepstrum vector corresponds to the input dgv vector
% this part is already functionized as ConvertData.m
% it can be used only when synthesis soon after getGMM_Qiao.m
% because loadGMM_Qiao.m doesn't define model.A
   Y2t = ConvertData(B', Model);
 sgram = cep2sgram(Y2t');


% for C programming, however, this part is re-written based on GMM_Mapping.m
%           B: input dgv data(interpolated), SNS * total_time
%    yMeanArr: \mu^B mean of target vector in mapping function (GMM), mixNo x dY    
%    xMeanArr: \mu^A mean of source vector in mapping function (GMM), mixNo x dX
%     xCovArr: \Sigma^{AA} covariance of source vector in mapping function, mixNo x dX x dY
% TinvxCovArr: \Sigma^{BA}\Sigma^{AA}^{-1}, mixNo x dX x dY
%        gNum: GMM mixture number
%        WArr: weights of GMM, 1 x mixNo
%

% dNum = size(B, 2); % the number of data
%  
% h = zeros(dNum, Model.gNum);
% for n = 1:dNum
%     for q = 1:Model.gNum
%         x = B(:, n)';
%         h1 = (x - Model.xMeanArr(q, :)) * inv(squeeze(Model.xCovArr(q, :, :))) * (x - Model.xMeanArr(q, :))';
%         h2 = (2*pi)^9 * det(squeeze(Model.xCovArr(q, :, :)))^(1/2);
%         h(n, q) = exp(-h1/2)/h2;
%     end
% end
% for n = 1:dNum
%     h(n, :) = h(n, :) .* Model.WArr;
% end
% sh = sum(h, 2);
% for n = 1:dNum
%     for q = 1:Model.gNum
%         h(n, q) = h(n, q) / sh(n, 1);
%     end
% end
% clear sh;
% 
% Y = zeros(dNum , Model.dY);
% for n = 1:dNum
%      y = zeros(1, Model.dY);
%      for q = 1:Model.gNum
%          y1 = squeeze(Model.TinvxCovArr(q, :, :)) * (B(:, n) - Model.xMeanArr(q, :)');
%          y2 = Model.yMeanArr(q, :) + y1';
%          y = y + y2 * h(n, q);
%      end
%      Y(n, :) = y;
% end
% 
% sgram = cep2sgram(Y');
 
 
%% resynthesis speech
% define sampling rate
fs = 16000; % 10[ms] = 100[Hz]

% make ap matrix
ap2 = sgram;

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
wavwrite(sy/32768, fs, 16, fout);