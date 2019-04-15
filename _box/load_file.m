% 2008-11-15
% load files
%

clear all

cd 'C:\Documents and Settings\Aki\My Documents'
A = importdata('avg_var.csv');

fmax = size(A.data, 1);

L = zeros(fmax/3, 2);
S = zeros(fmax/3, 2);
SL = zeros(fmax/3, 2);

for ii = 1:fmax
    if rem(ii, 3)==1
        L((ii+2)/3, 1) = A.data(ii, 1);
        L((ii+2)/3, 2) = A.data(ii, 2);
    elseif rem(ii, 3)==2
        S((ii+1)/3, 1) = A.data(ii, 1);
        S((ii+1)/3, 2) = A.data(ii, 2);
    elseif rem(ii, 3)==0
        SL(ii/3, 1) = A.data(ii, 1);
        SL(ii/3, 2) = A.data(ii, 2);
    end
end

fout_l = fopen('l.csv', 'wt');
fout_s = fopen('s.csv', 'wt');
fout_sl = fopen('sl.csv', 'wt');

row = fmax/15;
fprintf(fout_l, ' ,aiueo,euoai,ioaue,oaeiu,ueioa\n');
fprintf(fout_s, ' ,aiueo,euoai,ioaue,oaeiu,ueioa\n');
fprintf(fout_sl, ' ,aiueo,euoai,ioaue,oaeiu,ueioa\n');

for jj = 1:row;
    fprintf(fout_l, 'No. %d, %f, %f, %f, %f, %f\n', jj-1, L(jj, 1), L(jj+row, 1), L(jj+2*row, 1), L(jj+3*row, 1), L(jj+4*row));
    fprintf(fout_s, 'No. %d, %f, %f, %f, %f, %f\n', jj-1, S(jj, 1), S(jj+row, 1), S(jj+2*row, 1), S(jj+3*row, 1), S(jj+4*row));
    fprintf(fout_sl, 'No. %d, %f, %f, %f, %f, %f\n', jj-1, SL(jj, 1), SL(jj+row, 1), SL(jj+2*row, 1), SL(jj+3*row, 1), SL(jj+4*row));
end

%% release file pointer
fclose(fout_l);
fclose(fout_s);
fclose(fout_sl);