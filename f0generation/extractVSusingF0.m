%
% 2011/06/10
% extractVSusingF0.m extracts Voicing Strength using F0 information
% extracted using STRAIGHT
%
% INPUT
% finWav: input wav file, recorded in 16[kHz]
% finF0: input F0 file, extracted by finWav
% mode: 0: no F0 file, 1: use f0 file
%
% LINKS
% STRAIGHT, sgram2cep.m, cep2sgram.m
%
% NOTES
% this function is not checked
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%

clear all, fclose all, clc;

% test
finF0  = 'J:\!speech\ATR503\mht\f0\a01.f0';
finWav = 'J:\!speech\ATR503\mht\wav\mht_a_001.wav';
mode   = 1; 


%% definition
% f0: f0 countour, 1 x fnum
% F0 search Lower Bound : 40
f0LowerBound = 40;
% F0 search Upper Bound : 800
f0UpperBound = 800;
% F0 default Window Length: 40 (= 640[point])
f0WindowLength = 640;
% F0 frame Update Interval:  1 (=  16[point])
f0FrameUpdateInterval = 16;


%% load wav and f0 files
%  s: signal recorded with 16[kHz], bit num x 1
% fs: 16000 [Hz]
[s, fs] = wavread(finWav);

% f0: f0 countour, 1 x fnum
if mode == 0;
    [f0, ap, prmF0] = exstraightsource(s, fs);
else mode == 1;
    f0 = load(finF0);
    f0 = f0';
end

% for check
%[n3sgram, prmP] = exstraightspec(s, f0, fs);


%% calculate correlation
M       = length(f0); % the number of frames
wstart  = 0;
s       = [s; s(1:f0WindowLength)]; % consider the same signal will come after the signal
VS      = [];
for ii = 1:M
    wstart = (ii-1) * f0FrameUpdateInterval + 1;
    wsignal = s(wstart:(wstart + f0WindowLength - 1));
    if f0(ii) ~= 0
    % voiced region
        wlambda = floor(fs/f0(ii)); % wave length
        
        % get correlation for voicing strength
        w1    = s(wstart:wstart+wlambda-1);
        w2    = s(wstart+wlambda:wstart+2*wlambda-1);
        VSmat = corrcoef(w1, w2)+1; % make it non-negative
        VS    = [VS, VSmat(1, 2)];
    else
    % unvoiced region
        VS_ = 0;
        %BUF = [];
        for jj = f0LowerBound:f0UpperBound;
            wlambda = floor(fs/jj); % wave length
        
            % get correlation for voicing strength
            w1    = s(wstart:wstart+wlambda-1);
            w2    = s(wstart+wlambda:wstart+2*wlambda-1);
            VSmat = corrcoef(w1, w2)+1; % make it non-negative
            buf   = VSmat(1, 2);
            %BUF = [BUF, buf];
            if buf > VS_;
                VS_ = buf;
            end
        end
        VS = [VS, VS_*0.8]; % weight 0.8?
        %break;
    end
end

hold on
    plot(f0);
    plot(VS * 50, 'r');
hold off