function mapping(fin, fout, Model)
%
% 2009-10-04
% space mapping based on GMM
%
% INPUT
% dr: the directory where GMM data located
% fin: a input dgv file from 18 sensor CyberGlove
% fout: where wav file will be located
%
% LINK
% STRAIGHT, cep2sgram
%
% REMARKS
% Hand Motions are interpolated by every 1ms after loading
% f0 = 140[Hz] constant
% ap = -20 constant
% NUM = 22; % number of signal from dataglove
% SNS = 18; % number of sensor
% DEG = 17; % degree of cepstrum, usually SNS - 1
% JNT = 36; % degree per a joint vector
%
% REMARKS
% GMM data is made by HTK or getGMM.m
%



%% definition
NUM = 22;
SNS = 18;
DEG = 17;
JNT = 36;
 
% dX = 18; % dimensions of source vectors
% dY = 18; % dimensions of target vectors
% gNum = 8; % GMM mixture number
 
%dir = 'L:\hmts\091008_SP\model4-n\matrix_3-8_2';
%fin = 'L:\hmts\091008_SP\model4-n\dgv\open\_a1.dgv';
%fout = 'L:\hmts\091008_SP\model4-n\syn_3-4\_a1.wav';


%% get dgv vector from input 
finM = fopen(fin, 'rb');
 
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

fclose(finM);

Y2t = ConvertData(B', Model);


%% resynthesis speech
% get spectrogram from cepstrum time series
sgram = cep2sgram(Y2t');

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


%% manual

% %% load GMM data
% 
% fn_WArr = [dir '\WArr.txt'];
% fn_xMeanArr = [dir '\xMeanArr.txt'];
% fn_yMeanArr = [dir '\yMeanArr.txt'];
% fn_xCovArr = [dir '\xCovArr.txt'];
% fn_TinvxCovArr = [dir '\TinvxCovArr.txt'];
% 
% WArr = importdata(fn_WArr);
% xMeanArr_ = importdata(fn_xMeanArr);
% yMeanArr_ = importdata(fn_yMeanArr);
% xCovArr_ = importdata(fn_xCovArr);
% TinvxCovArr_ = importdata(fn_TinvxCovArr);
% 
% xMeanArr = zeros(gNum, dX);
% yMeanArr = zeros(gNum, dY);
% xCovArr = zeros(gNum, dX, dY);
% TinvxCovArr = zeros(gNum, dX, dY);
% 
% for r = 1:gNum
%     for c = 1:dX
%         xMeanArr(r, c) = xMeanArr_(gNum*(c-1)+r);
%     end
% end
% for r = 1:gNum
%     for c = 1:dY
%         yMeanArr(r, c) = yMeanArr_(gNum*(c-1)+r);
%     end
% end
% for r = 1:dX
%     for c = 1:dY
%         for g = 1:gNum;
%             xCovArr(g, c, r) = xCovArr_(gNum*dY*(r-1) + gNum*(c-1) + g);
%         end
%     end
% end
% for r = 1:dX
%     for c = 1:dY
%         for g = 1:gNum;
%             TinvxCovArr(g, c, r) = TinvxCovArr_(gNum*dY*(r-1) + gNum*(c-1) + g);
%         end
%     end
% end
% 
% 
% %% calculate f(a)
% %           B: input dgv data(interpolated), SNS * total_time
% %    yMeanArr: \mu^B mean of target vector in mapping function (GMM), mixNo x dY    
% %    xMeanArr: \mu^A mean of source vector in mapping function (GMM), mixNo x dX
% %     xCovArr: \Sigma^{AA} covariance of source vector in mapping function, mixNo x dX x dY
% % TinvxCovArr: \Sigma^{BA}\Sigma^{AA}^{-1}, mixNo x dX x dY
% %        gNum: GMM mixture number
% %        WArr: weights of GMM, 1 x mixNo
% 
% t = 1;
% % input
% x_ = B(:, t) / 1000;
% x = x_';
% 
% % Gaussian distribution with mean \mu^A, covariance \Sigma^{AA}     
% % N(x: xMeanArr_i, xCovArr_i)
% Gaus = zeros(gNum, 1);
% 
% for q = 1:gNum
%     % load the maen and the covariance
%     xMeanArr_q = xMeanArr(q, :);
%     xCovArr_q = [];
%     for r = 1:dX
%         xCovArr_q = [xCovArr_q; xCovArr(q, :, r)];
%     end
%     xCovArr_qi = inv(xCovArr_q);
%     % calculate Gaussian Distribution
%     Gaus(q, 1) = exp((x - xMeanArr_q)* inv(xCovArr_qi) * (x - xMeanArr_q)' * (-1/2)) / (((2 * pi) ^ (dX/2)) * sqrt(det(xCovArr_qi)));
% end
% 
% for q = 1:gNum
%     h(q, 1) = WArr(q, 1) * Gaus(q, 1) / dot(WArr, Gaus);
% end
% 
% 
% %% calculate f(a)
% %num = 1;
% 
% f = zeros(DEG + 1, 1);
% %     while num < mix_num + 1
% 
% % 変数の読み込み
% % 		mean_v = mean(num, 1:20);
% %         mean_w = mean(num, 21:40);
% %         covar_vv = covar(40 * (num - 1) + 1 : 40 * (num - 1) + 20, 1:20);
% %         covar_wv = covar(40 * (num - 1) + 21: 40 * num,            1:20);
%  
% %         h = alpha(num, 1) * gd(num, 1) / dot(alpha, gd);
% for q = 1:gNum
%     %f = f + h * (mean_w' + covar_wv * inv(covar_vv) * (v - mean_v)');
%     TinvxCovArr_q = [];
%     for r = 1:dX
%         TinvxCovArr_q = [TinvxCovArr_q; xCovArr(q, :, r)];
%     end
%     %TinvxCovArr_qi = inv(TinvxCovArr_q);
%     %f_ = yMeanArr(q, :)' + TinvxCovArr_q * (x - xMeanArr(q, :))';
%     %f = f + h(q, 1) * f_;
% end
% % 
% %         num = num + 1;
% %     end
% % 
% %     % ファイルへの書き出し
% %     fwrite(fod, A(1, t), 'float');
% %     fwrite(fod, f, 'float');
% %     
% % 	t = t + 1;
% % end
% % 
% % fclose(fid);
% % fclose(fod);