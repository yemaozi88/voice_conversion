% 2008-2-18
% check rank
%

clear all
cd /Users/Aki/voice_conversion/warping/5vowels/fftcep20deg/check_rank;
dirlist = dir('original/*.fftcep');
dirlength = length(dirlist);

%- deg of cepstrum
deg = 20;

if 0
%%%%% warping %%%%%
i = 1;
while i < dirlength + 1
	num = strrep(dirlist(i).name,'.fftcep','');
	fid = fopen(['original/' num '.fftcep'], 'rb');
	fod1 = fopen(['converted/warp_' num '.fftcep'], 'wb');
	fod2 = fopen(['joint_vector/joint_' num '.fftcep'], 'wb');
	fod3 = fopen(['upto100/joint_' num '.fftcep'], 'wb');
	A = fread(fid, [deg + 1, inf], 'float');

	%- number of flames
	fmax = length(A(1, :));

	%- warping matrix
	B = mk_matrix2(deg, 0.2);

	t = 1;
	while t < fmax + 1
		C = B * A(2:(deg + 1), t);
	
		%% write down warping data to file
		fwrite(fod1, A(1, t), 'float');
		fwrite(fod1, C, 'float');

		fwrite(fod2, A(2:(deg + 1), t), 'float');
		fwrite(fod2, C, 'float');
	
		%% write down joint vectors to file
		%% header for HTK
		%- number of samples in file
		fwrite(fod3, fmax, 'int');
		%- sample period in 100ms
		fwrite(fod3, 50000, 'int');
		%- number of bytes per sample
		fwrite(fod3, deg * 8, 'short');
		%- a code indicating the sample kind
		fwrite(fod3, 9, 'short');
				
		t = t + 1;
	end

	fclose(fid);
	fclose(fod1);
	fclose(fod2);
	fclose(fod3);
	i = i + 1;
end
%%%%%
end


%%%%% evaluation %%%%%
%% joint_vector.fftcep is combined all joint_*.fftcep
fid = fopen(['joint_vector/joint_vector.fftcep'], 'rb');
Z = fread(fid, [2 * deg, inf], 'float');

	%- number of flames
	n = length(Z(1, :));
	
	%- mean of Z
	mean_Z = mean(Z')';
	
	%- deviation of Z
	dev_Z = zeros(2 * deg, n);
	t = 1;
	while t < n + 1
		dev_Z(:, t) = Z(:, t) - mean_Z;
		t = t + 1;
	end
	
	%- covariance
	covar_Z = zeros(2 * deg, 2 * deg);
	i = 1;
	while i < 2 * deg + 1
		j = 1;
		while j < 2 * deg + 1
			covar_Z(i, j) = dev_Z(i, :) * dev_Z(j, :)' / n;
			j = j + 1;
		end
		i = i + 1;
	end
	
fclose(fid);
%%%%%