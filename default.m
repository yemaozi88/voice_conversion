%
% 2019/05/11
% default settings for the voice conversion project.
%
% AUTHOR
% Aki Kunikoshi
% a.kunikoshi@gmail.com
%

%% definition
clear all, fclose all, clc;

dirVC  = 'c:\OneDrive\Research\McRoberts\voice_conversion';
dirIn  = [dirVC '\fws'];
dirOut = [dirVC '\mht'];

% directory to be made.
dirNameList = {'wav', 'scep', 'resyn'};
dirKindList = {'In', 'Out'};
for dirName_ = dirNameList
    dirName = string(dirName_);
    for dirKind_ = dirKindList
        dirKind = string(dirKind_);

        folderName = sprintf('dir%s_%s', dirKind, dirName);
        evalStr = sprintf(...
            '%s = [dir%s ''\\%s''];', ...
            folderName, dirKind, dirName);
        eval(evalStr);
        
        evalStr = sprintf(...
            'isExist = exist(''%s'', ''dir'');', ...
            folderName);
        eval(evalStr);

        % if not exist, make the directory.
        if ~isExist
            mkdir(folderName);
        end % if
    end % dirKind
end % dirName
clear dirKind dirKind_ dirKindList
clear dirName dirName_ dirNameList
clear evalStr folderName isExist