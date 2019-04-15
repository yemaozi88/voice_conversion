%
% 2012/01/16
% do processes to all files in the directory
% 
% HISTORY
% this code is a simple version of doDir.m
%
% AUTHOR
% Aki Kunikoshi (D3)
% yemaozi88@gmail.com
%

fclose all, clear all, clc

%% definition
dirIn = 'C:\research\!speech\ATR503\mht\wav';

%% directory processing
dirlist = dir([dirIn '\*.*']);
dirlength = length(dirlist);

for ii = 1:dirlength
%for ii = 1:3 % test data is the first one
	% except ".", "..", "DS_Store"
  	if length(dirlist(ii).name) > 3
        filename = dirlist(ii).name;
        %filename = strrep(dirlist(ii).name, '.wav', '');
        %[pathstr, name, ext] = fileparts(ffeature);

        %% do
        %extractScep(fin, fout, 18);
    end
end % ii