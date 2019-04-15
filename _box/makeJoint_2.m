%%
%% 2008-2-14
%% make jointvector without DP matching
%%

clear all
cd /Users/Aki/voice_conversion/warping/5vowels/fftcep20deg/check_rank;
dirlist = dir('original');
dirlength = length(dirlist);

%- deg of cepstrum
deg = 20;

i = 1;
while i < dirlength + 1
	num = strrep(dirlist(i).name,'.fftcep','')

	%% except ".", "..", "DS_Store"
	if length(num) == 3

		fid1 = fopen(['original/' num '.fftcep'], 'rb');
		fid2 = fopen(['converted/warp_' num '.fftcep'], 'rb');
		fod = fopen(['gmm/upto100/joint_' num '.fftcep'], 'wb');

		A1 = fread(fid1, [deg + 1, inf], 'float');
		A2 = fread(fid2, [deg + 1, inf], 'float');

		%- number of flames
		fmax = length(A1(1, :));

	%% header for HTK
		%- number of samples in file
		fwrite(fod, fmax, 'int');
		%- sample period in 100ms
		fwrite(fod, 50000, 'int');
		%- number of bytes per sample
		fwrite(fod, deg * 8, 'short');
		%- a code indicating the sample kind
		fwrite(fod, 9, 'short');
	
		t = 1;
		while t < fmax + 1
	
		%% write to file
			fwrite(fod, A1(2:(deg + 1), t), 'float');
			fwrite(fod, A2(2:(deg + 1), t), 'float');
		t = t + 1;
		end

		fclose(fid1);
		fclose(fid2);
		fclose(fod);

	end
	i = i + 1;
end