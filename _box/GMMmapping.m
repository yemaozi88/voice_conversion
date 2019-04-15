%
% 2010-07-15
% GMM mapping for LSP + F0 + VS + their delta
%
% LINK
% loadBin.m, gen3weightmatrix_lspf0vsddd.m
%
% REMARKS
% This program is based on train24_2.m written by v-wedi
%
% AUTHOR
% Kunikoshi Aki, The University of Tokyo, Japan
% ( Microsoft Research Asia intern since May until August, 2010 )
%   Email: kunikoshi@gavo.t.u-tokyo.ac.jp
%  

clear all, clc;


%% definition

% parameters
vectorsize = 108;
sourcesize = 54;
lsporder = 24;
Sentence = 100;

% for models
x_idx=[1:52];
y_idx=[53:104];
s_idx=[1:26];
d_idx=[27:52];

% the location of directories
VCdir  = 'D:\users\v-akkuni\VC';

%% English M2F
%srcSpk = 'bdl';
%tarSpk = 'clb';
%trainRange = 101:200;
%testRange  = 201:220;

%% Chinese M2F
% srcSpk = 'LY';
% tarSpk = 'MS';
% trainRange = 1:100;
% testRange  = 101:120;

%% Chinese F2F
srcSpk = 'HF';
tarSpk = 'ZT';
trainRange = 100501:100600;
testRange  = 100701:100720;


% the directory name setup
srcVCdir = [VCdir '\' srcSpk];
tarVCdir = [VCdir '\' tarSpk];
convVCdir = [VCdir '\' srcSpk '2' tarSpk];


srcFeaDir = [VCdir '\' srcSpk '\lspf0vsddd'];
srcIF0Dir = [srcVCdir '\if0\'];

srcVSdir = [srcVCdir '\vs_smoothed'];
tarVSdir = [tarVCdir '\vs_smoothed'];

convOutDir = [convVCdir '\conv'];

tic;
for MixNum = [1, 2, 4, 8, 16, 32, 64, 128, 256]
MixNum = num2str(MixNum)


    % load GMM result, gparam2
    load(['lspf0vsddd_mix' MixNum]);
    % load srcMean, srcVar, tarMean and tarVar to calculate target pitch
    GaussDist_f0;


    for fNum = testRange;
        fNum = num2str(fNum);
        disp(fNum)

        srcFeaFile = [srcFeaDir '\' fNum '.lspf0vsddd'];
        srcIF0File = [srcIF0Dir '\' fNum '.if0'];


        %% F0 conversion using conventional method
        Pitch = load(srcIF0File);
        idx = find(Pitch>0);
        Pitch(idx) = (tarVar/srcVar)*(Pitch(idx) - srcMean) + tarMean;
        clear idx;


        %% Spectral comversion

        %load files
        %disp(srcFeaFile)
        if  ~exist(srcFeaFile, 'file')
            disp(sprintf('Error: file not exists : %d\n', fNum));
            continue;
        end

        srcFea = loadBin(srcFeaFile, 'float', sourcesize);
        Frame_Num = size(srcFea, 2);    
        gain = srcFea(25,:);
        srcFea = srcFea([1:24, 26:27, 28:51, 53:54], :);

        % adjust source form to GMM training data, like dtwProc4lspf0vsddd.m as follows:
        %  1-24: LSP x 1
        %   (25: gain)
        %    26: f0 x 0.1
        %    27: Voicing Strength x 1
        % 28-51: delta LSP x 50
        %   (52: delta gain)
        %    53: delta f0 x 10
        %    54: delta VS x 10
        weight = [ones(1, lsporder), 0.1, 1, ones(1, lsporder)*50, 10, 10]';
        weight = repmat(weight, 1, Frame_Num);
        srcFea = srcFea .* weight;
        clear weight;

        % DIM = sourcesize
        % N = Frame_Num
        [DIM, N] = size(srcFea);

        % extract source part of gparam2
        parm_src = gmix_strip(gparam2, x_idx);

        % Version of lqr_evp for transposed input data (DIM-by-nsamp)
        %   it returns log mode PDF in columns, that is the output of lqr_eval only, 
        %   not including mixing weights for each column.
        lpyk = lqr_evp(parm_src, srcFea, 1); 

        % get the likelyhood and to find the max component
        %  choose the most likely distribution
        %  instead of taking weighted summation over all distribution
        srcmax = max(lpyk, [], 2);

        for i=1:N,
            k = find(lpyk(i,:)==srcmax(i));
            tmpvar = gparam2.modes(k).cholesky_covar;  % Dim *Dim Dim = vectorsize=100;
            tmpvar = tmpvar' * tmpvar; % covariance matrix

            sxx = 1./diag(tmpvar(x_idx, x_idx)); 
            % in fact it is the inversed Sxx 
            %  = diagonal part of inv(diagonal matrix of tmpvar(x_idx, x_idx))
            syy = diag(tmpvar(y_idx, y_idx));
            sxy = diag(tmpvar(x_idx, y_idx));

            mean = gparam2.modes(k).mean(y_idx) ...
                   + sxy.*sxx.*(srcFea(:,i) - gparam2.modes(k).mean(x_idx));
            Sigma = syy - sxy.*sxx.*sxy;

            miuM_1(:,i) = mean(s_idx); % mean of spectral part
            miuM_2(:,i) = mean(d_idx); % mean of delta spectral part
            Sigma_1(:,i) = Sigma(s_idx); % sigma of spectral part
            Sigma_2(:,i) = Sigma(d_idx); % sigma of delta spectral part
        end
        clear sxx sxy syy;
        clear tmpvar mean;
        clear parm_src srcFea;


        % creat Em Dm and then W matirx, 
        % weight factor here are with doubt, i think.
        % DIM = 52;
        DIM = DIM/2;
        % DIM = 26, for static and delta feature both
        miuVec_1 = reshape(miuM_1, DIM*N, 1);
        miuVec_2 = reshape(miuM_2, DIM*N, 1);

        SigmaVec_1 = reshape(Sigma_1, DIM*N, 1);
        SigmaVec_2 = reshape(Sigma_2, DIM*N, 1);

        clear miuM_1 miuM_2 Sigma;
        clear Sigma_1 Sigma_2;

        % in gene2weightmatrix_2, we donnot denerate the gain volunm.
        gen3weightmatrix_lspf0vsddd; % gen3weightmatrix_2 with comments

        SigmaM_1 = spalloc(DIM*N, DIM*N, DIM*N);
        SigmaM_2 = spalloc(DIM*N, DIM*N, DIM*N);

        for Mdim = 1 : LPCORDER*N
             SigmaM_1(Mdim, Mdim) = SigmaVec_1(Mdim);
             SigmaM_2(Mdim, Mdim) = SigmaVec_2(Mdim);
        end

        invSigmaM_1 = inv(SigmaM_1);
        invSigmaM_2 = inv(SigmaM_2);
        clear SigmaVec_1 SigmaVec_2 ;

    %     W_2 = 100*W_2;
        traW_invSigmaM_1 = W_1' * invSigmaM_1;
        traW_invSigmaM_2 = W_2' * invSigmaM_2;

        A_1 = traW_invSigmaM_1*W_1;
        A_2 = traW_invSigmaM_2*W_2;

        b_1 = traW_invSigmaM_1*miuVec_1;
        b_2 = traW_invSigmaM_2*miuVec_2;

        tarFea = (A_1+A_2)\(b_1+b_2); 
        clear A_1 A_2 b_1 b_2 traW_invSigmaM_1 traW_invSigmaM_2 invSigmaM_1 invSigmaM_2 ;
        clear miuVec_1 miuVec_2;
        %DIM = 24
        tarFea = reshape(tarFea, DIM, length(tarFea)/DIM);

        % mediem filter
        for vec=1:DIM
            tarFea(vec,:) = medfilt1(tarFea(vec,:));
        end
        tarFea(25, :) = medfilt1(tarFea(25, :), 7); % F0
        Pitch = medfilt1(Pitch);

        TAR = zeros(DIM + 1, N);
        for ii=2:N+1,
            %resort lsp feature
            tarFea(1:24, ii-1) = sort(tarFea(1:24, ii));
            TAR(1:24, ii-1) = tarFea(1:24, ii); % LSP
            TAR(25, ii-1) = gain(ii-1);         % gain
            TAR(26, ii-1) = tarFea(25, ii);     % F0
            TAR(27, ii-1) = tarFea(26, ii);     % VS
        end

        % calculate F0
    %     tarF0 = zeros(1, N);
    %     for ii = 1:N
    %         if TAR(27, ii) > 0.8
    %             tarF0(ii) = exp(TAR(26, ii));
    %         end
    %     end


        %% Output   
        outfeafile   = [convOutDir '\' fNum '_mix' MixNum '.lsp'];
        outf0LCfile  = [convOutDir '\' fNum '_mix' MixNum '.f0lc'];
        outf0Rawfile = [convOutDir '\' fNum '_mix' MixNum '.f0raw'];
        %outf0VSfile  = [convOutDir '\' fNum '_mix' MixNum '.f0vs'];
        outVSfile    = [convOutDir '\' fNum '_mix' MixNum '.vsraw'];

        fp_fea   = fopen(outfeafile, 'w');
        fp_f0lc  = fopen(outf0LCfile, 'wt');
        fp_f0raw = fopen(outf0Rawfile, 'wt');
        %fp_f0vs  = fopen(outf0VSfile, 'wt');
        fp_vs    = fopen(outVSfile, 'wt');

        for ii = 1:N
            fwrite(fp_fea, TAR(1:25, ii), 'float32');
            fprintf(fp_f0lc, '%f\n', Pitch(ii));
            fprintf(fp_f0raw, '%f\n', exp(TAR(26, ii)));
            %fprintf(fp_f0vs, '%f\n', tarF0(ii));
            fprintf(fp_vs, '%f\n', TAR(27, ii));
        end

        fclose(fp_fea);
        fclose(fp_f0lc);
        fclose(fp_f0raw);
        %fclose(fp_f0vs);
        fclose(fp_vs);

        %clear tarFea;
    end

    %% convert VS considering Global Variance(GV)
    disp('convert VS considering GV');
    convVS;

end % end for MixNum
toc;