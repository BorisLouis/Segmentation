clear
close all
clc

idx2File = 1;
pathForTest = 'N:\Project Z-stacks\New stacks\1K_0.25mg_30 nm_sample 1\tif_1K_0.25 mg-30 nm_Sample 1\SegmentedStacks';
path = uigetdir;

%%
imSegmentation.check(path, idx2File);
