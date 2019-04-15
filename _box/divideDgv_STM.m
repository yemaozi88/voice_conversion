%
% 2009-10-03 
% check dgv transition
%
% LINK
% loadDgv.m
%


%% definition
NUM = 22; % number of signal from dataglove
SNS = 18; % number of sensor
DEG = 18; % number of the deg of scep, coef deg + 1 (e.g. 0-18 -> DEG = 19),

excB = 220;
excE_ = 220;

L = 2;


%% write out
%dir = 'O:\hmts\100430_IS\n16\dgv\n';
%nn = 2;
%cv = 'nu';


%% load dgv
%A_ = loadDgv('K:\hmts\091008_SP\model4-n\dgv\n\3\no.dgv');
%filename = [dir '\' num2str(nn) '\' cv '.dgv'];
filename = 'O:\hmts\100430_IS\n16\dgv\open\no3.dgv';
A_ = loadBin(filename, 'uchar', 26);
fmax_ = size(A_, 2);


%% forming 
% the first dgv data is wrong
A = zeros(SNS, fmax_ - 1);
%A(1:SNS, :) = A_(5:SNS+4, :);
A(1:SNS, :) = A_(5:SNS+4, 2:fmax_);
fmax = size(A, 2);

excE = fmax - excE_;


%% get the difference
% Y = zeros(SNS, fmax-1);
% for ii=1:fmax-1
%     Y(:, ii) = X(:, ii+1) - X(:, ii);
% end


%% calculate a_i(m) = C/B
% B = \sum^L_{n=-L}n^2
B = 0;
for n = -L:1:L
    B = B + n^2;
end

% C = \sum^L_{n=-L}MFCC_i(n+m)
STM = zeros(1, fmax);
for m = 1:L
    STM(m) = 0;
end
for m = fmax-L+1:fmax
    STM(m) = 0;
end

for m = 1+L:fmax-L
    E = 0;
    for dim = 1:DEG-1
        A_ = A(dim + 1, :);
        C = 0;
        for n = -L:1:L
            C = C + A_(n+m) * n;
        end
        E = E + (C/B)^2;
    end
    STM(m) = E;
end

plot(STM)


%% detect max STM point (= the transition between phones)
% first
[STMmaxB, STMmaxBt_] = max(STM(1:excB));
STMmaxBt = STMmaxBt_;

% middle
%[STMmaxM, STMmaxMt_] = max(STM(excB+1:fmax-excE));
%STMmaxMt = STMmaxMt_ + excB;

% end
[STMmaxE, STMmaxEt_] = max(STM(fmax-excE+1:fmax));
STMmaxEt = STMmaxEt_ + fmax - excE;

disp(STMmaxBt)
%disp(STMmaxMt)
disp(STMmaxEt)


%% devide scep file
%Ac = A(:, STMmaxt-1);
%Av = A(:, STMmaxt:fmax);


%% output 2 scep files
%manDivideDgv(filename, STMmaxBt, STMmaxEt);
% filenamec = [dirOut '\' num2str(nn) '\' cv '-c.dgv'];
% filenamev = [dirOut '\' num2str(nn) '\' cv '-v.dgv'];
% 
% foutc = fopen(filenamec, 'wb');
% foutv = fopen(filenamev, 'wb');
% 
% fwrite(foutc, Ac, 'uchar');
% fwrite(foutv, Av, 'uchar');
% 
% fclose(foutc);
% fclose(foutv);