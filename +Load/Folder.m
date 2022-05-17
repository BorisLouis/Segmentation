% The aim of this code is to allow the user to choose a folder to analyze
% while precising which type of files he wants to analyze in the folder.
% The function also create a directory to store output and output the path
% to this folder. 

%The function return the list of files of the folder that correspond to 
%the encoded extension and the path to the current folder and the output 
%folder.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boris Louis (https://github.com/BorisLouis)                             %
% Website : Boris Louis: https://borislouis.github.io/                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [file2Analyze,currentFolderName,outDir] = Folder(fileExt,output)

%Open UI to select folder
mainFolderName = uigetdir;
assert(ischar(mainFolderName),'User canceled the selection of file, excecution aborted');


%extract the name of the current folder
idx = strfind(mainFolderName,filesep) ;
currentFolderName = mainFolderName(idx(end)+1:end) ;

%Remove dots from the name
currentFolderName = regexprep(currentFolderName,'\.','_');

% generate folder to store output
outDir = [mainFolderName filesep output];
status = mkdir(outDir);

%Extract the part of the folder that is a tif file
Folder_Content = dir(mainFolderName);
index2Images   = contains({Folder_Content.name},fileExt);
file2Analyze = Folder_Content(index2Images);

if isempty(file2Analyze)
    warning('No %s file found in the selected directory')
end

end