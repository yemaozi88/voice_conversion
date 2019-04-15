% 2009-50-2
% do processes to all files in the directory
% 


%% directory processing
% 16 gestures which i can make easily
%ges = [4];
%ges = [4, 9, 13, 14, 15, 16, 21, 22, 25, 27];
%ges = [13, 14, 15, 16, 21, 22, 25, 27];
gmax = size(ges, 2);

jj = 1;
while  jj < gmax + 1
    
    numI = ges(1, jj);    % integer
    numS = num2str(numI); % string
    if length(numS) == 1  % if numS is 1 digit then add 0 at the head
        numS = ['0' numS];
    end

    %% definition
    % PCA4joint
    dirIn = ['C:\hmts\100430_INTERSPEECH\n' numS '\joint\joint-cmv'];
    dirOut = ['C:\hmts\100430_INTERSPEECH\n' numS '\joint\joint-cmv_PCA-8'];

    % reduceData
    dirOut2 = ['C:\hmts\100430_INTERSPEECH\n' numS '\joint\joint-cmv_PCA-8_reduced-10'];

    % getGMM
    dirOut3_1 = ['C:\hmts\100430_INTERSPEECH\n' numS '\matrix_1'];
    dirOut3_2 = ['C:\hmts\100430_INTERSPEECH\n' numS '\matrix_2'];
    dirOut3_4 = ['C:\hmts\100430_INTERSPEECH\n' numS '\matrix_4'];
    dirOut3_8 = ['C:\hmts\100430_INTERSPEECH\n' numS '\matrix_8'];

    % mappingPCA2GMM_Qiao
    dirIn4 = ['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\test'];
    dirOut4_1 = ['C:\hmts\100430_INTERSPEECH\n' numS '\syn_1'];
    dirOut4_2 = ['C:\hmts\100430_INTERSPEECH\n' numS '\syn_2'];
    dirOut4_4 = ['C:\hmts\100430_INTERSPEECH\n' numS '\syn_4'];
    dirOut4_8 = ['C:\hmts\100430_INTERSPEECH\n' numS '\syn_8'];


%Model = loadGMM_Qiao('C:\hmts\100419_ProgressReport\n16\matrix-cmv_2_2', 8, 18, 2, 5);
dirEval = ['C:\hmts\100430_INTERSPEECH\n' numS '\Eval'];
[Evec, Eval, u] = loadEval(dirEval);


%dirlist = dir(num2str(nn));
dirlist = dir(dirIn);
dirlength = length(dirlist);

ii = 1;
while ii < dirlength + 1
	% except ".", "..", "DS_Store"
  	if length(dirlist(ii).name) > 3 
        filename = dirlist(ii).name;
        
        % PCA4joint
        fin = [dirIn '\' filename];
        fout = [dirOut '\' filename];
        %disp(['PCA4joint:' fin])
        %PCA4joint(fin, fout, Evec, u, 8);
        
        % reduceData
        fout2 = [dirOut2 '\' filename];        
        %disp(['reduceData:' fout2])
        %reduceData(fout, fout2, 'float', 26, 10);
    end
    ii = ii + 1;
end
clear dirlist;

    % synthesis speech
    getGMM_Qiao(dirOut2, dirOut3_1, 1);
    Model = loadGMM_Qiao(dirOut3_1, 8, 18, 1, 5);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\1\na.dgv'], [dirOut4_1 '\closed_na1.wav'], Model, Evec, u, 8);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\2\na.dgv'], [dirOut4_1 '\closed_na2.wav'], Model, Evec, u, 8);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\3\na.dgv'], [dirOut4_1 '\closed_na3.wav'], Model, Evec, u, 8);

    getGMM_Qiao(dirOut2, dirOut3_2, 2);
    Model = loadGMM_Qiao(dirOut3_2, 8, 18, 2, 5);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\1\na.dgv'], [dirOut4_2 '\closed_na1.wav'], Model, Evec, u, 8);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\2\na.dgv'], [dirOut4_2 '\closed_na2.wav'], Model, Evec, u, 8);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\3\na.dgv'], [dirOut4_2 '\closed_na3.wav'], Model, Evec, u, 8);

    getGMM_Qiao(dirOut2, dirOut3_4, 4);
    Model = loadGMM_Qiao(dirOut3_4, 8, 18, 4, 5);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\1\na.dgv'], [dirOut4_4 '\closed_na1.wav'], Model, Evec, u, 8);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\2\na.dgv'], [dirOut4_4 '\closed_na2.wav'], Model, Evec, u, 8);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\3\na.dgv'], [dirOut4_4 '\closed_na3.wav'], Model, Evec, u, 8);
    
    getGMM_Qiao(dirOut2, dirOut3_8, 8);
    Model = loadGMM_Qiao(dirOut3_8, 8, 18, 8, 5);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\1\na.dgv'], [dirOut4_8 '\closed_na1.wav'], Model, Evec, u, 8);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\2\na.dgv'], [dirOut4_8 '\closed_na2.wav'], Model, Evec, u, 8);
    mappingPCA2GMM_Qiao(['C:\hmts\100430_INTERSPEECH\n' numS '\dgv\n\3\na.dgv'], [dirOut4_8 '\closed_na3.wav'], Model, Evec, u, 8);

    jj = jj + 1;
end

%fclose(fout);