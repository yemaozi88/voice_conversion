clear all

load_gmm

%- 混合数 mix_num
%- 重み alpha (mix_num, 1)
%- 平均 mean (mix_num, 40)
%- 共分散行列 covar (40 * mix_num, 40)

%fid = fopen('../../../../Users/Aki/voice_conversion/mht2fws/original/mht/m_f_001.mcep', 'r');
%fod = fopen('../../../../Users/Aki/voice_conversion/mht2fws/converted/gmm.mcep', 'wb');
fid = fopen('../../../../Users/Aki/voice_conversion/mht2fws_2/original/mht/m_a_001.mcep', 'r');
fod = fopen('../../../../Users/Aki/voice_conversion/mht2fws_2/converted/gmm8_a_001.mcep', 'wb');


A = fread(fid, [21, inf], 'float');

%- 最後のフレーム数
max = length(A(1, :));

t = 1;
while t < max + 1
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

    % ファイルへの書き出し
    fwrite(fod, A(1, t), 'float');
    fwrite(fod, f, 'float');
    
	t = t + 1;
end

fclose(fid);
fclose(fod);