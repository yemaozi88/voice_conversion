%
% 2008/10/17 
% Principal Component Analysis of dgv data from CyberGlove
%
% HISTORY
% 2010-04-21 cut out getPCA part and renamed PCA_dgv into exPCA
% 2010-04-20 functionized loadBinDir 
%

clear all;

NUM = 22; % number of signal from dataglove
SNS = 18; % number of sensor


%% calculate eigen vectors & eigen values for PCA
%getPCA;


%% load eigen vectors & eigen values for PCA
[Evec, Eval, u] = loadEval('C:\hmts\100419_ProgressReport\Eval_all');


%% contribution ratio
% Eval_total = 0;
% for i = 1:SNS;
%     Eval_total = Eval_total + Eval(i, 1);
% end
% 
% % accumulative summation and its potion
% Eval_accum = 0;
% cr = zeros(SNS, 1);
% for i = 1:SNS;
%     Eval_accum = Eval_accum + Eval(i, 1);
%     cr_ = Eval_accum / Eval_total * 100;
%     cr(i, 1) = cr_;
% end

%plot(cr, 'ko-');


%% 28gestures
dirIn = 'C:\hmts\100419_ProgressReport\ex3';

% 16 gestures which i can make easily
ges = [1, 2, 4, 8, 9, 11, 13, 14, 15, 16, 21, 22, 25, 27, 28];
gmax = size(ges, 2);

hold on

ii = 1;
while  ii < gmax + 1
    numI = ges(1, ii);     % integer
    numS = num2str(numI); % string
    if length(numS) == 1  % if numS is 1 digit then add 0 at the head
        numS = ['0' numS];
    end
    filename1 = [dirIn '\' numS '-' numS '_1.dgv'];
    filename2 = [dirIn '\' numS '-' numS '_2.dgv'];
    filename3 = [dirIn '\' numS '-' numS '_3.dgv'];
    X1 = loadBin(filename1, 'uchar', 26);
    X2 = loadBin(filename2, 'uchar', 26);
    X3 = loadBin(filename3, 'uchar', 26);
    
    X = [X1, X2, X3];
    
    clear X1;
    clear X2;
    clear X3;

    X = X(5:22, :); % extract 18 sensor data
    Xm = mean(X');
    Y = PCA_Trans(Xm, Evec, u, 2);
    
    %disp(numS);
    %disp(Y);
    %visuallize the location of 28gestures
    if strcmp(numS, '13') == 1
        plot(Y(:,1), Y(:,2), 'b+');
    else
        plot(Y(:,1), Y(:,2), 'b+');
    end

    ii = ii + 1;
end


   
%% transition
tr = '14-02';
tr = num2str(tr);

ftr1 = ['C:\hmts\100430_INTERSPEECH\dgv-all\' tr '_1.dgv'];
%ftr2 = ['C:\hmts\100430_INTERSPEECH\dgv-all\' tr '_2.dgv'];
%ftr3 = ['C:\hmts\100430_INTERSPEECH\dgv-all\' tr '_3.dgv'];

%ftr1 = ['C:\hmts\PCA\dgv\16gestures\dgv\' tr '_1.dgv'];
%ftr2 = ['C:\hmts\PCA\dgv\16gestures\dgv\' tr '_2.dgv'];
%ftr3 = ['C:\hmts\PCA\dgv\16gestures\dgv\' tr '_3.dgv'];

X1 = loadBin(ftr1, 'uchar', 26);
X1_fmax = size(X1, 2);
X1 = X1(5:22, 2:X1_fmax)'; % extract 18 sensor data
Y1 = PCA_Trans(X1, Evec, u, 2);

% X2 = loadBin(ftr2, 'uchar', 26);
% X2_fmax = size(X2, 2);
% X2 = X2(5:22, 2:X2_fmax)'; % extract 18 sensor data
% Y2 = PCA_Trans(X2, Evec, u, 2);
% 
% X3 = loadBin(ftr3, 'uchar', 26);
% X3_fmax = size(X3, 2);
% X3 = X3(5:22, 2:X3_fmax)'; % extract 18 sensor data
% Y3 = PCA_Trans(X3, Evec, u, 2);

plot(Y1(:,1), Y1(:,2), 'r.');
%plot(Y2(:,1), Y2(:,2), 'g.');
%plot(Y3(:,1), Y3(:,2), 'b.');

hold off