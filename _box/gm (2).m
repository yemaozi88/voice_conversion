% 2008-4-28
% make gaussian model
%

clear all
[x, fs] = wavread('/Users/Aki/voice_conversion/gy2mym/gm/g_1_01.wav');
[f0raw, ap, prmF0] = exstraightsource(x, fs);
fs2 = fs/1.8;

%- deg of cepstrum
deg = 20;

%% load joint.vector
fid_a = fopen('/Users/Aki/voice_conversion/gy2mym/joint/joint.scep', 'rb');
A = fread(fid_a, [2*deg, inf], 'float');

m = mean(A')';

C = cov(A');
iC = inv(C);


%% load input dgv2 file
fid_b = fopen('/Users/Aki/voice_conversion/gy2mym/gm/g_1_01.scep','rb');
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
    scep(1, t) = B(1, t);
    scep(2:deg+1, t) = f;

    t = t + 1;
end

sgram = cep2sgram(scep);
    
%% write to wav file
sy = exstraightsynth(f0raw, sgram, ap, fs2);
wavwrite(sy/32768, fs2, 16, '/Users/Aki/voice_conversion/gy2mym/gm/gm_1_01.wav');

fclose(fid_a);
fclose(fid_b);
%fclose(fod);