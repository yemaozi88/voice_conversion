
%
% 2009-08-10 
% Mix speech of a consonant and a vowel
%
% LINK
% STRAIGHT
%
% REMARKS
% wav files should be recorded 16000[Hz]
%


%% definition
finC = 'C:\hmts\091008_SP\concatenate\parts\na-c.wav';
finV = 'C:\hmts\091008_SP\concatenate\parts\aa_l.wav';
fout = 'C:\hmts\091008_SP\concatenate\syn\na-c_na-c_aa_l_30-60-60_2.wav';
damp_ms = 30;
enhance_ms = 60;
overlap_ms = 60;
consonant_rt = 1;
damp_rt = 0.9;
enhance_rt = 0.9;


%% load wav files
[c, fs_c] = wavread(finC);
[v, fs_v] = wavread(finV);
if fs_c ~= fs_v
    error('Sampling frequency of 2 files are not the same.');
end

% data size
size_c = size(c, 1);
size_v = size(v, 1);

disp(['the length of the consonant: ' num2str(size_c) '[bit] (' num2str(size_c/16) '[ms])'])
disp(['the length of the vowel: ' num2str(size_v) '[bit] (' num2str(size_v/16) '[ms])'])

% damping part of the end of a consonant
damp_bit = damp_ms * 16;
if damp_bit > size_c
    error('overlapped part is longer than the consonant part.');
end

% enhancing part of the top of a vowel
enhance_bit = enhance_ms * 16;
if enhance_bit > size_v
    error('overlapped part is longer than the vowel part.');
end

% overlapped part
overlap_bit = overlap_ms * 16;

disp(['the damping part: ' num2str(damp_bit) '[bit] (' num2str(damp_ms) '[ms])'])
disp(['the enhancing part: ' num2str(enhance_bit) '[bit] (' num2str(enhance_ms) '[ms])'])
disp(['the overlapped part: ' num2str(overlap_bit) '[bit] (' num2str(overlap_ms) '[ms])'])


%% windowing (the blackman window to keep continuity)
% damp the end of a consonant
x_c = 0.5:1/(2*(damp_bit - 1)):1;
y_c = 0.42 - 0.5 * cos(2 * pi * x_c) + 0.08 * cos(4 * pi * x_c);
for i = (size_c - damp_bit + 1):size_c
     c(i, 1) = y_c(1, i - (size_c - damp_bit)) * c(i, 1);
end
% for i = 1:size_c
%     c(i, 1) = 0.8 * c(i, 1);
% end

% enhance the top of a vowel
x_v = 0:1/(2*(enhance_bit - 1)):0.5;
y_v = 0.42 - 0.5 * cos(2 * pi * x_v) + 0.08 * cos(4 * pi * x_v);
for i = 1:enhance_bit
      v(i, 1) = y_v(1, i) * v(i, 1);
end


%% combine the consonant part and the vowel part
mix = zeros(size_c + size_v - overlap_bit, 1);

% consonant part
for i = 1:(size_c - overlap_bit)
     mix(i, 1) = c(i, 1) * consonant_rt;
end
 
% overlapped part 
for i = (size_c - overlap_bit + 1):size_c
     mix(i, 1) = c(i, 1) * damp_rt + v(i - (size_c - overlap_bit), 1) * enhance_rt;
end

% vowel part
for i = (size_c + 1):(size_c + size_v - overlap_bit)
     mix(i, 1) = v(i - (size_c - overlap_bit), 1);
end
 
 
%% write to wav file
% for check, use soundsc(mix, fs_c)
wavwrite(c, fs_v, 16, fout);