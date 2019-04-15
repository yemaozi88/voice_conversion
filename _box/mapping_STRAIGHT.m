clear all

load_gmm

%- 混合数 mix_num
%- 重み alpha (mix_num, 1)
%- 平均 mean (mix_num, 40)
%- 共分散行列 covar (40 * mix_num, 40)

[x, fs] = wavread('../../../../Users/Aki/voice_conversion/STRAIGHT-cep/wav/mht/m_a_001.wav');
[f0raw, ap, prmF0] = exstraightsource(x, fs);

[n3sgram, prmP] = exstraightspec(x, f0raw, fs);

fid = fopen('../../../../Users/Aki/voice_conversion/STRAIGHT-cep/scep/big_endian/mht/m_a_001.scep-b', 'r');

A = fread(fid, [21, inf], 'float');


%- 最後のフレーム数
frame = length(A(1, :));

%- 
B = zeros(21, frame);

t = 1;
while t < frame + 1
    %- 入力
    v = A(2:21, t)';

    %% 平均mean_v, 分散covar_vvのガウス分布（Gaussian Distribution）を計算する
    
    % i番目の混合におけるガウス分布を格納した行列
    gd = zeros(mix_num, 1);

    num = 1;
    while num < mix_num + 1
        % 変数の読み込み
        mean_v = mean(num, 1:20);
        covar_vv = covar(40 * (num - 1) + 1:40 * (num - 1) + 20, 1:20);

        % ガウス分布の計算
        gd(num, 1) = exp((v - mean_v)*inv(covar_vv) * (v - mean_v)' * (-1/2)) / (((2 * pi) ^ 10) * sqrt(det(covar_vv)));
        num = num + 1;
    end

    %% F(v)を計算する
    num = 1;
    f = zeros(20, 1);
    while num < mix_num + 1
        % 変数の読み込み
		mean_v = mean(num, 1:20);
        mean_w = mean(num, 21:40);
        covar_vv = covar(40 * (num - 1) + 1 : 40 * (num - 1) + 20, 1:20);
        covar_wv = covar(40 * (num - 1) + 21: 40 * num,            1:20);

        % F(v)の計算
        h = alpha(num, 1) * gd(num, 1) / dot(alpha, gd);
        f = f + h * (mean_w' + covar_wv * inv(covar_vv) * (v - mean_v)');

        num = num + 1;
    end
	
	B(1, t) = A(1, t);
	B(2:21, t) = f;
    %fwrite(fod, A(1, t), 'float');
    %fwrite(fod, f, 'float');
    
	t = t + 1;
end


sgram = cep2sgram(A);
f0raw2 = f0raw * 1.98;

%%%%% synthesis sound %%%%%
[sy, prmS] = exstraightsynth(f0raw, sgram, ap, fs);

%%%%% write to wav file %%%%%
wavwrite(sy/32768, fs, 16, '../../../../Users/Aki/voice_conversion/STRAIGHT-cep/converted/m_a_001.gmm4.wav');

fclose(fid);