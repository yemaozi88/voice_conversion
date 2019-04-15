%
% 2009/10/03 
% getSTMpeak.m checks dgv transition
%
% LINK
% loadBin.m, getSTM.m
%
% HISTORY
% 2011/04/05 linked to getSTM.m
% 
% AUTHOR
% Aki Kunikoshi (D1)
% yemaozi88@gmail.com
%


%% definition
fstart  = 1;  % file starts
excB    = 150; % a consonant ends
excE_   = 150; % a vowel starts
fend    = 250; % file ends

%fin = 'J:\ProbabilisticIntegrationModel\distortion\bestModel\1113\conventionalModel\dgv\1\08-01.dgvs';
fin = '/Volumes/RESEARCH/ProbabilisticIntegrationModel/distortion/bestModel/1113/conventionalModel/n/3/no.dgvs';
type = 'dgvs';


%% output
%dir = 'O:\hmts\100430_IS\n16\dgv\n';
% dir = 'J:\ProbabilisticIntegrationModel\distortion\bestModel\1113\conventionalModel\n-cmv';
% nn = 2;
% cv = 'nu';
% type = '.dgv';


%% get STM
STM = getSTM(fin, type);
fmax = size(STM, 2);
excE = fmax - excE_;

plot(STM)


%% detect max STM point (= the transition between phones)
% first
[STMmaxB, STMmaxBt_] = max(STM(fstart:excB));
STMmaxBt = STMmaxBt_ + fstart;

% middle
%[STMmaxM, STMmaxMt_] = max(STM(excB+1:fmax-excE));
%STMmaxMt = STMmaxMt_ + excB;

% end
if fmax < fend
    fend = fmax;
end
[STMmaxE, STMmaxEt_] = max(STM(fmax-excE+1:fend));
STMmaxEt = STMmaxEt_ + fmax - excE;

disp(['consonant ends at ' num2str(STMmaxBt) '[frame]'])
%disp(STMmaxMt)
disp(['vowel starts at ' num2str(STMmaxEt) '[frame]'])


%% devide scep file
manDivide(fin, STMmaxBt, STMmaxEt, 'dgvs');