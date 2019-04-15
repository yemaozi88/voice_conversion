% 2008-4-28
% make gaussian model
%

clear all
[x, fs] = wavread('C:\Users\kunikoshi\Documents\HMTS\080424arai_model\5vowel\cutout\wav\001.wav');
[f0raw, ap, prmF0] = exstraightsource(x, fs);

%- deg of cepstrum
deg = 22;

%% load joint.vector
fid_a = fopen('C:\Users\kunikoshi\Documents\HMTS\080424arai_model\5vowel\joint\5vowel.joint', 'rb');
A = fread(fid_a, [2*deg, inf], 'float');

m = mean(A')';

C = cov(A');
iC = inv(C);


%% load input dgv file
fid_b = fopen('C:\Users\kunikoshi\Documents\HMTS\080424arai_model\5vowel\cutout\dgv\001.dgv','rb');
B = fread(fid_b, [deg + 1, inf], 'float');
max = length(B(1,:));
scep = zeros(deg + 1, max);

%fod = fopen('/Users/Aki/voice_conversion/HMTS-mac/kuni_model/gm/a-i.scep','wb');

t = 1;
while t < max + 1
%    f = m(deg+1:2*deg) + C(deg+1:2*deg, 1:deg)/C(1:deg, 1:deg) * (B(:, t) - m(1:deg));
%    scep(:, t) = f;
%    fwrite(fod, f, 'float');

    f = m(deg+1:2*deg) + C(deg+1:2*deg, 1:deg)/C(1:deg, 1:deg) * (B(2:deg+1, t) - m(1:deg));
    scep(1, t) = A(1, 100);
    scep(2:deg+1, t) = f;

    t = t + 1;
end

sgram = cep2sgram(scep);
    
%% write to wav file
sy = exstraightsynth(f0raw, sgram, ap, fs);
wavwrite(sy/32768, fs, 16, 'C:\Users\kunikoshi\Documents\HMTS\080424arai_model\gm001.wav');

fclose(fid_a);
fclose(fid_b);
%fclose(fod);