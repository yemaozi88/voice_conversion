%
% 2011/09/28
% simple version of MLtest4MixNum.m
%
% LINKS
% loadBinDir.m, trainGMM.m
%
% HISTORY
% 2011/10/25 cleaned up the code for Keigo
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%

clc, clear all, fclose('all');

%% definition
% the directory for training data
dirIn_train = 'J:\H2Swith16deg_0243\joint\jointS2H_11_1of8';
% the directory for testing data
dirIn_test  = 'J:\H2Swith16deg_0243\joint\jointS2H_22_1of8';
% output gmm object and the log file
dirOut      = 'J:\H2Swith16deg_0243\joint\S2Hmodel_withDelta_full_withPCA';
% maximum mixture number
gMax_       = 1024;

% joint
Y = loadBinDir(dirIn_train, 'float', 32); % load all binary data in dirIn_train
U = loadBinDir(dirIn_test, 'float', 32);  % load all binary data in dirIn_test
Y = Y';
U = U';


%% preparation for ML calculation
gNumArray = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
gMax = length(gNumArray);
gMax = log2(gMax_)+1;
ML = zeros(gMax, 1);
clear gMax_

% log
flog = fopen([dirOut '\log.txt'], 'wt');
fprintf(flog, '<Maximum Likelihood>\n');
fprintf(flog, 'The number of data : %d\n\n', size(Y, 1));

for g = 1:gMax
tic
    gNum = gNumArray(g);
    disp(gNum)

    obj = trainGMM(Y, gNum, 0); % 0 - full, 1 - diagonal
    fOut = [dirOut '\S2Hmodel_mix' num2str(gNum) '_obj'];
    save(fOut, 'obj');
    %load(fOut)
    clear fOut
        
    L_ = pdf(obj, U);
    L = log(max(L_));
    ML(g, 1) = L;
    disp([num2str(gNum) ':' num2str(L)]);

    % log
    fprintf(flog, 'The number of mixtures : %d\n', gNum);
    %fprintf(flog, '\ttrain : %6.4f\n', obj.NlogL);
    fprintf(flog, '\ttest  : %6.4f\n', L);
    
    save([dirOut '\ML'], 'ML');
toc
end
fclose(flog);
clear flog