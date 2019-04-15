%
% 2011/08/09
% VoiceConversion.m estimates target speaker's speech from input speaker's
% speach based on space mapping
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%


%% voice conversion
clear all, fclose all, clc;
TYPE = 1; % static only: 0, with delta:1


%% definition
% f0 information
dirSrcF0_train = 'J:\VoiceConversion\STRAIGHT-based\fws\f0_train';

dirSrcF0train = 'J:\VoiceConversion\STRAIGHT-based\fws\f0_train';
dirTgtF0train = 'J:\VoiceConversion\STRAIGHT-based\mht\f0_train';

dirSrcF0 = 'J:\VoiceConversion\STRAIGHT-based\fws\f0_test';
dirTgtF0 = 'J:\VoiceConversion\STRAIGHT-based\mht\f0_test';  

if TYPE == 0
    % static only
    dirScep = 'J:\VoiceConversion\STRAIGHT-based\fws\scep18_test';
    dirOut  = 'J:\VoiceConversion\STRAIGHT-based\fws2mht_reduced1of5\converted';

    % gmm information
    load('J:\VoiceConversion\STRAIGHT-based\fws2mht_reduced1of5\gmm_full\jointModel_mix128_obj');
elseif TYPE == 1
    % with delta
    dirScep = 'J:\VoiceConversion\STRAIGHT-based\fws\scep18_test';
    dirFea  = 'J:\VoiceConversion\STRAIGHT-based\fws\scep18+delta_test';
    dirOut  = 'J:\VoiceConversion\STRAIGHT-based\fws2mht_withDelta_reduced1of5\converted';

    % gmm information
    load('J:\VoiceConversion\STRAIGHT-based\fws2mht_withDelta_reduced1of5\gmm_full\jointModel_mix128_obj');
end


%% directory processing
for n = 1:10
tic
    if n < 10
        nStr = ['0' num2str(n)];
    else
        nStr = num2str(n);
    end
    nStr = ['j' nStr];
    disp(nStr)


%% load files
    if ismac == 1
        finSrcF0 = [dirSrcF0 '/' nStr '.f0'];
        finTgtF0 = [dirTgtF0 '/' nStr '.f0'];
        if TYPE == 1
            % with delta
            finFea = [dirFea '/' nStr '.feature'];
        end
        finScep  = [dirScep '/' nStr '.scep'];
        foutScep = [dirOut '/' nStr '.scep'];
        foutWav  = [dirOut '/' nStr '.wav'];
    else
        finSrcF0 = [dirSrcF0 '\' nStr '.f0'];
        finTgtF0 = [dirTgtF0 '\' nStr '.f0'];
        if TYPE == 1
            % with delta
            finFea = [dirFea '\' nStr '.feature'];
        end
        finScep  = [dirScep '\' nStr '.scep'];
        foutScep = [dirOut '\' nStr '.scep'];
        foutWav  = [dirOut '\' nStr '.wav'];
    end
    srcF0 = load(finSrcF0);
    tgtF0 = load(finTgtF0);
        
    srcF0train = loadDir(dirSrcF0train);
    tgtF0train = loadDir(dirTgtF0train);

    srcScep  = loadBin(finScep, 'float', 19);
    if TYPE == 0
        % static only
        srcScep_ = srcScep(2:19, :); % remove energy
    elseif TYPE == 1
        % with delta
        srcFea  = loadBin(finFea, 'float', 36);
    end


%% F0 conversion
    fnum = length(srcF0);
    for ii = 1:fnum
        srcF0_2(ii, 1) = (srcF0(ii, 1) - mean(srcF0train)) * sqrt(var(tgtF0train)/var(srcF0train)) + mean(tgtF0train);
    end
    clear srcF0 tgtF0 srcF0train tgtF0train


%% spectral conversion
    if TYPE == 0
        % static only
        %tgtScep = gmmvc(srcScep_, obj);
        tgtScep = gmmMapping(srcScep_, obj);
        % using energy of source speaker
        tgtScep = [srcScep(1, :); tgtScep];
    elseif TYPE == 1
        % with delta
        tgtScep = gmmMappingWithDelta2(srcFea, obj);
        % modify the data length
        %tgtScep(:, fnum) = [];
        %tgtScep(:, 1) = [];
        
        % Saito's code
%         it = 1;
%         tgtScep = gmmvc_delta(srcFea, obj, it);
%         tgtScep = tgtScep{1, it};
        
        % using energy of source speaker
        tgtScep = [srcScep(1, :); tgtScep(1:18, :)];
    end



%% write to scep file
    fScep = fopen(foutScep, 'wb');
    fNum = size(tgtScep, 2);
    for ii = 1:fNum
        fwrite(fScep, tgtScep(:, ii), 'float');
    end
    fclose(fScep);


%% synthesis speech using STRAIGHT
    % parameters for STRAIGHT synthesis
    sgram   = cep2sgram(tgtScep);
    ap      = repmat(-20, 513, fnum);
    fs      = 16000;

    % synthesis sound
    [sy, prmS] = exstraightsynth(srcF0_2, sgram, ap, fs);

    % write to wav file
    wavwrite(sy/32768, fs, 16, foutWav);
toc
end %n