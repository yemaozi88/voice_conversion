% 2008-4-28
% extract feature vector of Hand Motion from joint vector
%

%- degree of feature vector (= 1/2 degree of joint vector)
deg = 22;

fid = fopen('/Users/Aki/voice_conversion/HMTS-mac/kuni_model/joint-b/a-i.joint-b', 'rb');
fod = fopen('/Users/Aki/voice_conversion/HMTS-mac/kuni_model/gm/a-i.dgv2', 'wb');

A = fread(fid, [2 * deg, inf], 'float');

max = length(Z(1,:));

t = 1;
while t < max + 1
    fwrite(fod, A(1:22, t), 'float');    
    t = t + 1;
end

fclose(fid);
fclose(fod);