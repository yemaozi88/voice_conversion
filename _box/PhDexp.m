%
% HITORY
% 2011/07/03
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%
%

% signal
%x1 = -pi:0.01:pi;
%y1 = sin(x1);
%x2 = -1/2*pi:0.01:3/2*pi;
%   y2 = sin(x2);

clear all, fclose all, clc

%dirIn = '/Users/kunikoshi/Desktop/PhDexp';
dirIn = 'C:\Documents and Settings\kunikoshi\My Documents\Dropbox\document\PhDexp';


%% load data

if ismac == 1
    sonic1m = load([dirIn '/raw/1m.log']);
    sonic2m = load([dirIn '/raw/2m.log']);
    sonic3m = load([dirIn '/raw/3m.log']);
else
    sonic1m = load([dirIn '\raw\1m.log']);
    sonic2m = load([dirIn '\raw\2m.log']);
    sonic3m = load([dirIn '\raw\3m.log']);
end

for deg = 10:10:50
    if ismac == 1
        fdegP = [dirIn '/raw/+' num2str(deg) 'deg.log'];
        fdegM = [dirIn '/raw/-' num2str(deg) 'deg.log'];
    else
        fdegP = [dirIn '\raw\+' num2str(deg) 'deg.log'];
        fdegM = [dirIn '\raw\-' num2str(deg) 'deg.log'];
    end

    % plus
    if exist(fdegP)
        tmpP = sprintf('sonicP%ddeg=load(fdegP);', deg);
        eval(tmpP);
    else
        disp([fdegP 'does not exsit.']);
    end
    
    % minus
    if exist(fdegM)
        tmpM = sprintf('sonicM%ddeg=load(fdegM);', deg);
        eval(tmpM);
    else
        disp([fdegM 'does not exsit.']);
    end
end
clear fdegP fdegM tmpP tmpM deg

% second 0 deg
if ismac == 0
    sonic0deg  = load([dirIn '/raw/0deg.log']);
    sonic0deg2 = load([dirIn '/raw/0_2deg.log']);
else
    sonic0deg  = load([dirIn '\raw\0deg.log']);
    sonic0deg2 = load([dirIn '\raw\0_2deg.log']);
end
clear dirIn


%% extract appropriate period and re-scale from micron to mm

sonic1m = sonic1m(55:1054)/1000;    % 1854 x 1
sonic2m = sonic2m(26:1025)/1000;   % 1025 x 1
sonic3m = sonic3m(1001:2000)/1000;  % 2009 x 1

sonic0deg  = sonic0deg(7:1006)/1000;
sonic0deg2 = sonic0deg2(27:1026)/1000;

sonicP10deg = sonicP10deg(117:1116)/1000;
sonicP20deg = sonicP20deg(7:1006)/1000; 
sonicP30deg = sonicP30deg(9:1008)/1000; 
sonicP40deg = sonicP40deg(12:1011)/1000; 
sonicP50deg = sonicP50deg(7:1006)/1000; 

sonicM10deg = sonicM10deg(18:1017)/1000; 
sonicM20deg = sonicM20deg(7:1006)/1000; 
sonicM30deg = sonicM30deg(18:1017)/1000; 
sonicM40deg = sonicM40deg(7:1006)/1000; 
sonicM50deg = sonicM50deg(9:1008)/1000; 


%% get means
% distance
disMean = zeros(3, 1);
disMean(1) = mean(sonic1m);
disMean(2) = mean(sonic2m);
disMean(3) = mean(sonic3m);

% angle
angMean = zeros(12,1);
angMean(1)  = mean(sonicM50deg);
angMean(2)  = mean(sonicM40deg);
angMean(3)  = mean(sonicM30deg);
angMean(4)  = mean(sonicM20deg);
angMean(5)  = mean(sonicM10deg);
angMean(6)  = mean(sonic0deg);
angMean(7)  = mean(sonicP10deg);
angMean(8)  = mean(sonicP20deg);
angMean(9)  = mean(sonicP30deg);
angMean(10) = mean(sonicP40deg);
angMean(11) = mean(sonicP50deg);
angMean(12) = mean(sonic0deg2);


%% draw graphs
dis1m = sonic1m - disMean(1);
dis2m = sonic2m - disMean(2);
dis3m = sonic3m - disMean(3);

angM50deg = sonicM50deg - angMean(1);
angM40deg = sonicM40deg - angMean(2);
angM30deg = sonicM30deg - angMean(3);
angM20deg = sonicM20deg - angMean(4);
angM10deg = sonicM10deg - angMean(5);
ang0deg   = sonic0deg - angMean(6);
angP10deg = sonicP10deg - angMean(7);
angP20deg = sonicP20deg - angMean(8);
angP30deg = sonicP30deg - angMean(9);
angP40deg = sonicP40deg - angMean(10);
angP50deg = sonicP50deg - angMean(11);
ang0deg2 = sonic0deg2 - angMean(12);

%hist(ang0deg2, 20)


%% modify data
dis3m_ = [];
angM30deg_ = [];
angM20deg_ = [];
angP10deg_ = [];
for ii = 1:1000
    if dis3m(ii) < 0.1
        dis3m_ = [dis3m_; dis3m(ii)]; 
    end
    
    if angM30deg(ii) > -0.15
        angM30deg_ = [angM30deg_; angM30deg(ii)];
    end

    if angM20deg(ii) < 0.15
        angM20deg_ = [angM20deg_; angM20deg(ii)];
    end
    
    if angP10deg(ii) > 0 % or angP10deg(ii) < 10
        angP10deg_ = [angP10deg_; angP10deg(ii)];
    end    
end
%hist(angM30deg_, 20)
clear ii

disMean(3) = mean(dis3m_) + disMean(3);
angMean(3) = mean(angM30deg_) + angMean(3);
angMean(7) = mean(angP10deg_) + angMean(7);
%plot(angMean(1:11)/1000)

disVar = zeros(3, 1);
disVar(1)  = var(dis1m);
disVar(2)  = var(dis2m);
disVar(3)  = var(dis3m_);

angVar = zeros(12,1);
angVar(1) = var(angM50deg);
angVar(2) = var(angM40deg);
angVar(3) = var(angM30deg_);
angVar(4) = var(angM20deg_);
angVar(5) = var(angM10deg);
angVar(6) = var(ang0deg);
angVar(7) = var(angP10deg_);
angVar(8) = var(angP20deg);
angVar(9) = var(angP30deg);
angVar(10) = var(angP40deg);
angVar(11) = var(angP50deg);
angVar(12) = var(ang0deg2);
%plot(angVar(1:11)


%% real distance
% distance
x1   = disMean(2)-disMean(1);
dis1 = sqrt(4+0.04^2)-sqrt(1+0.04^2) - x1/1000;
x2   = disMean(3)-disMean(1);
dis2 = sqrt(9+0.04^2)-sqrt(1+0.04^2) - x2/1000;
clear x1 x2

% angle
a = atan(1/2);
TR = [];
for theta = -50:10:50;
    theta = theta * pi/180;
    TR_ = 1.0836 - 0.0416*sqrt(5)*cos(a+theta);
    TR_  = sqrt(TR_);
    TR = [TR; TR_];
end
clear TR_

angMean = angMean(1:11);
angVar  = angVar(1:11);

idealTR = zeros(11, 1);
realTR  = zeros(11, 1);
errorTR = zeros(11, 1);
for ii = 1:11
    idealTR(ii) = TR(ii) - TR(6);
    realTR(ii)  = (angMean(ii) - angMean(6))/1000;
    errorTR(ii) = idealTR(ii) - realTR(ii);
end
[idealTR, realTR, errorTR]