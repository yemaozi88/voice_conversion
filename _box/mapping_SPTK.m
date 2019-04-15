clear all

load_gmm

%- ������ mix_num
%- �d�� alpha (mix_num, 1)
%- ���� mean (mix_num, 40)
%- �����U�s�� covar (40 * mix_num, 40)

%fid = fopen('../../../../Users/Aki/voice_conversion/mht2fws/original/mht/m_f_001.mcep', 'r');
%fod = fopen('../../../../Users/Aki/voice_conversion/mht2fws/converted/gmm.mcep', 'wb');
fid = fopen('../../../../Users/Aki/voice_conversion/mht2fws_2/original/mht/m_a_001.mcep', 'r');
fod = fopen('../../../../Users/Aki/voice_conversion/mht2fws_2/converted/gmm8_a_001.mcep', 'wb');


A = fread(fid, [21, inf], 'float');

%- �Ō�̃t���[����
max = length(A(1, :));

t = 1;
while t < max + 1
    %- ����
    v = A(2:21, t)';

    %% ����mean_v, ���Ucovar_vv�̃K�E�X���z�iGaussian Distribution�j���v�Z����
    
    % i�Ԗڂ̍����ɂ�����K�E�X���z���i�[�����s��
    gd = zeros(mix_num, 1);

    num = 1;
    while num < mix_num + 1
        % �ϐ��̓ǂݍ���
        mean_v = mean(num, 1:20);
        covar_vv = covar(40 * (num - 1) + 1:40 * (num - 1) + 20, 1:20);

        % �K�E�X���z�̌v�Z
        gd(num, 1) = exp((v - mean_v)*inv(covar_vv) * (v - mean_v)' * (-1/2)) / (((2 * pi) ^ 10) * sqrt(det(covar_vv)));
        num = num + 1;
    end

    %% F(v)���v�Z����
    num = 1;
    f = zeros(20, 1);
    while num < mix_num + 1
        % �ϐ��̓ǂݍ���
		mean_v = mean(num, 1:20);
        mean_w = mean(num, 21:40);
        covar_vv = covar(40 * (num - 1) + 1 : 40 * (num - 1) + 20, 1:20);
        covar_wv = covar(40 * (num - 1) + 21: 40 * num,            1:20);

        % F(v)�̌v�Z
        h = alpha(num, 1) * gd(num, 1) / dot(alpha, gd);
        f = f + h * (mean_w' + covar_wv * inv(covar_vv) * (v - mean_v)');

        num = num + 1;
    end

    % �t�@�C���ւ̏����o��
    fwrite(fod, A(1, t), 'float');
    fwrite(fod, f, 'float');
    
	t = t + 1;
end

fclose(fid);
fclose(fod);