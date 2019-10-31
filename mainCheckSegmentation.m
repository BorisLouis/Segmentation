clear
close all
clc

idx2File = 1;%if several segmented data is are in the same folder
path = uigetdir;
Frame2Skip = 10; %to scroll around the stack faster

%%
imSegmentation.check(path, idx2File,10);
