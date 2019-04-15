% 2009-09-23
% divide using STM a scep file to a consonant part and a vowel part
%

clear all

%% definition
dirScep = 'L:\hmts\091008_SP\scep\n';
%dirScep = 'L:\scep';
dirOut = 'L:\hmts\091008_SP\scep\n-cv';
dirWav = 'L:\voice\consonant_suzuki\16k\n';
nn = 2;
cv = 'nu';
excB = 30;
excE = 50;

DEG = 19; % number of the deg of scep, coef deg + 1 (e.g. 0-18 -> DEG = 19),
L = 2; % the number of frames on each side of the current frame

% B = \sum^L_{n=-L}n^2
B = 0;
for n = -L:1:L
    B = B + n^2;
end


%% manual processing
filenameScep = [dirScep '\' num2str(nn) '\' cv '.scep'];
%filenameScep = 'L:\scep\test2.scep';
filenameWav = [dirWav '\' num2str(nn) '\' cv '.wav'];
%filenameWav = 'L:\wav\test2.wav';
        
fin = fopen(filenameScep, 'rb');

        
%% file read       
A = fread(fin, [DEG, inf], 'float');
fmax = size(A, 2);
fclose(fin);


%% calculate a_i(m) = C/B
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


%% detect max STM point (= the transition between phones)
% first
[STMmaxB, STMmaxBt_] = max(STM(1:excB));
STMmaxBt = STMmaxBt_;

% middle
[STMmaxM, STMmaxMt_] = max(STM(excB+1:fmax-excE));
STMmaxMt = STMmaxMt_ + excB;

% end
[STMmaxE, STMmaxEt_] = max(STM(fmax-excE+1:fmax));
STMmaxEt = STMmaxEt_ + fmax - excE;

% manual
% STMmaxBt = 25;
% STMmaxMt = 175;
% STMmaxEt = 355;
 
disp(STMmaxBt)
disp(STMmaxMt)
disp(STMmaxEt)


%% comparison between wavform and STM
[x, fs] = wavread(filenameWav);
bmax = size(x, 1); % maximam bit number

% adjust scale of x -> STM
wform = zeros(1, floor(bmax/16)); % convert bit -> msec
for jj = 1:size(x, 1)
   if rem(jj, 16) == 0
       wform(jj/16) = x(jj);
   end
end
hold on
plot(wform)
plot(STM*200,'LineWidth', 2, 'Color',[.6 0 0])
hold off

% adjust scale of STM -> x
% STM2 = zeros(bmax, 1);
% for jj = 1:bmax-16
%     STM2(jj) = STM(floor(jj/16)+1);
% end
% hold on
% plot(x)
% plot(STM2*200,'LineWidth', 2, 'Color',[.6 0 0])
% hold off


%% devide scep file
Ac = A(:, STMmaxBt:STMmaxMt-1);
Av = A(:, STMmaxMt:STMmaxEt);

 
%% output 2 scep files
filenamec = [dirOut '\' num2str(nn) '\' cv '-c.scep'];
filenamev = [dirOut '\' num2str(nn) '\' cv '-v.scep'];

foutc = fopen(filenamec, 'wb');
foutv = fopen(filenamev, 'wb');

fwrite(foutc, Ac, 'float');
fwrite(foutv, Av, 'float');

fclose(foutc);
fclose(foutv);