clear all

load_gmm

%- ������ mix_num
%- �d�� alpha (mix_num, 1)
%- ���� mean (mix_num, 40)
%- �����U�s�� covar (40 * mix_num, 40)

[x, fs] = wavread('../../../../Users/Aki/voice_conversion/STRAIGHT-cep/wav/mht/m_a_001.wav');
[f0raw, ap, prmF0] = exstraightsource(x, fs);

[n3sgram, prmP] = exstraightspec(x, f0raw, fs);

fid = fopen('../../../../Users/Aki/voice_conversion/STRAIGHT-cep/scep/big_endian/mht/m_a_001.scep-b', 'r');

A = fread(fid, [21, inf], 'float');


%- �Ō�̃t���[����
frame = length(A(1, :));

%- 
B = zeros(21, frame);

t = 1;
while t < frame + 1
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