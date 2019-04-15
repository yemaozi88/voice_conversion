% 2008-8-15
% get mean vector and covarience matrix of Gaussian Model

clear all


%% definition
NUM = 22; % number of signal from dataglove
SNS = 18; % number of sensor
DEG = 17; % degree of cepstrum, usually SNS - 1
JNT = SNS + DEG + 1; % degree per a joint vector

dirIn = 'L:\hmts\091008_SP\model4-n';
dirJoint = 'joint5';
dirOut = 'L:\hmts\091008_SP\model4-n\matrix_5-1';


%% load joint vector
str = ['cd ' dirIn];
eval(str);
dirlist = dir(dirJoint);
dirlength = length(dirlist);


%% count total length of joint vectors
tic
i = 1;
j = 1;
jlen = 0; % total length of joint vectors

while i < dirlength + 1
 	%% except ".", "..", "DS_Store"
   	if length(dirlist(i).name) > 3 
        filename = dirlist(i).name;
     	fin = fopen([dirIn '\' dirJoint '\' filename], 'rb');
        A_ = [fread(fin, [JNT, inf], 'float')];
        jlen = jlen + length(A_(1, :)) - 1; % the first col of every file is 0

        clear A_;
        fclose(fin);
        %disp(j)
        j = j + 1;
    end
    i = i + 1;
end
t1 = toc;


%% set joint vectors to A
tic
A = zeros(JNT, jlen);

i = 1;
j = 1;
f_start = 1; % (col number of A) - 1
while i < dirlength + 1
 	%% except ".", "..", "DS_Store"
    	if length(dirlist(i).name) > 3 
            filename = dirlist(i).name;
            fin = fopen([dirIn '\' dirJoint '\' filename], 'rb');
            A_ = [fread(fin, [JNT, inf], 'float')];
            
            fnum = length(A_(1, :));
            f_end = f_start + fnum - 2;
            
            %disp(f_start)
            %disp(f_end)

            A(:, f_start:f_end) = A_(:,2:fnum);

            f_start = f_start + fnum - 1;
            
            clear A_;
            fclose(fin);
            disp(j)
            j = j + 1;
     end
     i = i + 1;
end
t2 = toc;


%% normalize size of joint matrix, dgv part / 1000
%B = zeros(JNT, jlen);
%B(1:SNS, :) = A(1:SNS, :)/100;
%B(SNS+1:JNT, :) = A(SNS+1:JNT, :);


%% get mean vector and covariance matrix
tic
m = mean(A')';
C = cov(A');
t3 = toc;


%% write out matrix data
tic
fout_mv = fopen([dirOut '\mean_v.txt'], 'wt');
fout_mw = fopen([dirOut '\mean_w.txt'], 'wt');
fout_c = fopen([dirOut '\covar_mul.txt'], 'wt');

fprintf(fout_mv, '%f\n', m(1:SNS));
fprintf(fout_mw, '%f\n', m(SNS+1:JNT));
fprintf(fout_c, '%f\n', C(SNS+1:JNT, 1:SNS)/C(1:SNS, 1:SNS));

fclose(fout_mv);
fclose(fout_mw);
fclose(fout_c);

t4 = toc;