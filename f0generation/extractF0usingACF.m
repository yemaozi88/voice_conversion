%
% 2011/06/15
% extractF0usingACF.m extract F0 by calculating autocorrelation function (ACF)
%
% NOTES
% - For detailed information, refer following paper:
% "A Robust Algorithm for Pitch Tracking (RAPT)"
%   David Talkin, 1995 Elsevier Science B.V.
%
% HISTORY
% 2011/06/29 modified so that variables correspond to those in extractF0usingACF.m
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%

%clear all, fclose all, clc;

% test
%finWav  = 'J:\!speech\JapaneseConsonants\b\suzuki\16k\wav\1\bu.wav';
%finF0  = 'J:\ProbabilisticIntegrationModel\F0generation\speech\b\f0\1\bu.f0';
finWav = 'J:\!speech\ATR503\suzuki\16k_from_48k\wav\a01.wav';
finF0  = 'J:\!speech\ATR503\suzuki\16k_from_48k\f0\a01.f0';


%% load wav file
%  s: signal recorded with 16[kHz], bit num x 1
% fs: 16000 [Hz]
% f0: f0 countour, 1 x fnum
%   F0 search Lower Bound : 40
%   F0 search Upper Bound : 800
%   F0 default Window Length: 40 (= 640[point])
%   F0 frame Update Interval:  1 (=  16[point])
[s, fs] = wavread(finWav);
f0 = load(finF0);
clear finWav


%% definition
% T: sampling interval
T = 1/fs;

% ii: frame index (= i)
% F0 search Lower Bound : 40
f0LowerBound = 40;
% F0 search Upper Bound : 800
f0UpperBound = 800;

% t: analysis frame interval [s]
%   In STRAIGHT, F0 frame Update Interval: 1[ms] (= 16[point])
t = 10^(-3);
% z: analysis frame interval [samples]
z = t/T;

% w: analysis window size [s]
%   w can be chosen to be on the order of a single average glottal period
%   In STRAIGHT, F0 default Window Length: 40[ms] (= 640[point])
% idx = find(f0 > 0);
% fnumVoiced = length(idx);
% meanF0 = sum(f0)/fnumVoiced;
% w = 1/meanF0;
w = 40 * 10^(-3);
% n: analysis window size [samples]
n = w/T;
% M: sample number
M = round(length(s)/fs * 1000);
% K: lag, K<n
% 40-800 [Hz]
%K = 10;
KupperBound = fs/f0LowerBound;
KlowerBound = fs/f0UpperBound;

clear fs T t w
%clear f0LowerBound f0UpperBound


%% calculate Voising Strength at lag k and analysis frame i
%s = [s; s];
s_ = repmat(0, n, 1); % add one frame at the beginning
s = [s; s_];
clear s_

%ACF = zeros(K, M);
ACF = zeros(KupperBound-KlowerBound+1, M);
%for k = 0:K-1;
for k = KlowerBound:KupperBound  
    for ii = 0:M-1
        m = ii * z;
        
        ACF_ = 0;
        for jj = m:(m+n-k-1)
            ACF_ = ACF_ + s(jj+1) * s(jj+k+1);
        end
        ACF(k+1, ii+1) = ACF_;
    end
end


%% get F0
idx = [];
f02 = [];
for ii = 1:M;
    [ACF2_, idx_] = max(ACF(:, ii));
    idx   = [idx, idx_];
    if idx_ == KlowerBound+1
        f02_ = 0;
    else
        f02_ = 16000/idx_;
        if f02_ > F0upperBound
            f02_ = tmp;
        end
    end
    f02 = [f02, f02_];
    tmp = f02_;
end
%f02 = medfilt1(f02, 15);
clear idx idx_ tmp f02_ KupperBound KlowerBound

% shift f02 to the half length of the window ahead
f03 = f02;
gap = n/z; % n: analysis window size ?
gap = gap/2;
for ii = 1:M;
    if ii > gap;
        f03(ii) = f02(ii-gap);
    else
        %f03(ii) = f02(ii+M-gap);
        %f03(ii) = f02(ii);
        f03(ii) = 0;
    end
end
%f03 = medfilt1(f03, 10);
clear n z gap
clear f02


%% check
% check
s4plot = [];
for ii = 1:length(s)
    if rem(ii, 16) == 1
        s4plot = [s4plot, s(ii)];
    end
end
hold on
    plot(f0);
    plot(f03, 'r');
    plot(s4plot*2000, 'k');
hold off