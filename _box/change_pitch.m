% 2008-2-
% change pitch
%

clear all

cd /Users/Aki/voice_conversion/mht2fws/original;
dirlist = dir('mht/*.pitch');
dirlength = length(dirlist);

i = 1;
while i < dirlength + 1

	num = dirlist(i).name;

	fid = fopen(['mht/' num], 'rb');
	fod = fopen(['2.2' num], 'wb');

	A = fread(fid, [1, inf], 'float');

	%- number of flames
	max = length(A(1, :));

	t = 1;
	while t < max + 1

		f = A(1, t) / 2.2;

		% write out to file
		fwrite(fod, f, 'float');

		t = t + 1;
	end

	fclose(fid);
	fclose(fod);
	
	i = i + 1;
end