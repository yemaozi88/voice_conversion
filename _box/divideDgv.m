%function divideDgv(fin)
% divide a dgv file to a consonant part and a vowel part
%
% divideDgv(fin)
%
% INPUT
% fin: a input dgv file
% fout: write out for the same folder to fin, consonant part -> fin-c.wav, vowel part -> fin-v.wav
%
% LINK
%
% REMARKS
% NUM: the number of dgv degrees
%
% HISTORY
% 2009-09-21 functionized
% 

clear all

%% definition
DEG = 26; % number of the deg of dgv, time: int(4 byte) + 18 sensors data: uchar(1 byte) + 4 zeros
dir = 'K:\hmts\091008_SP\dgv';
dgvb = 'uu'; % dgv data used for a border
dgvc = 'b'; % dgv data for a consonant
dgvv = 'o'; % dgv data for a vowel

% i = 1;
% while i < dirlength + 1
% 	% except ".", "..", "DS_Store"
%   	if length(dirlist(i).name) > 3 
%         filename = strrep(dirlist(i).name,'.scep','');
%         disp(filename)
%     end
%     i = i + 1;
% end

%% load data of a border gesture
B = [];
for n = 1:3
    filename = [dir '\vowels\' num2str(n) '\' dgvb '.dgv'];
    str = ['fb = fopen(''' filename ''', ''rb'');'];
    eval(str)
    
    B_ = fread(fb, [DEG, inf], 'uchar');
    B = [B, B_(:, 2:size(B_, 2))];
end
fclose(fb);


%% get mean vector of a border gesture data
avg = mean(B')';


%% load input data
%for n = 1:3;
n = 3;
    filename = [dir '\' dgvc '\' num2str(n) '\' dgvc dgvv '.dgv'];
    fin = fopen(filename, 'rb');

    A_ = fread(fin, [DEG, inf], 'uchar');
    A = A_(:, 2:size(A_, 2));
    tmax = size(A, 2);


%% calculate the difference between input data and average
    for t = 1:tmax
         dif_ = A(:, t) - avg;
         dif(1, t) = dot(dif_, dif_);
    end
    difmin = min(dif);
    difmint = 0;
    plot(dif)

    % compare with dgvb and dgvv, same:1, different:0
    cmp = strcmp(dgvb, strcat(dgvv, dgvv));
    th = (max(dif) - min(dif)) / 30;
    for t = 1:tmax
        if cmp == 0
            if difmin == dif(1, t)
                difmint = t;
                break;
            end
        else
            if difmin + th >= dif(1, t)
                difmint = t;
                break;
            end
        end
    end
    disp(difmint)


%% devide dgv file
    Ac = A(:, 1:difmint);
    Av = A(:, difmint+1:tmax);


%% output 2 dgv files
    filenamec = [dir '\' dgvc '-cv\' num2str(n) '\' dgvc dgvv '-c.dgv'];
    filenamev = [dir '\' dgvc '-cv\' num2str(n) '\' dgvc dgvv '-v.dgv'];

    foutc = fopen(filenamec, 'wb');
    foutv = fopen(filenamev, 'wb');

    fwrite(foutc, Ac, 'uchar');
    fwrite(foutv, Av, 'uchar');

    fclose(foutc);
    fclose(foutv);
%end


%% check
% fp = fopen('K:\hmts\091008_SP\dgv\b-cv\1\bu-c.dgv', 'rb');
% C = fread(fp, [DEG, inf], 'uchar');