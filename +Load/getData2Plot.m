% The aim of this function is to extract from the data output of
% mainPoreSizeCalc some specific properties (prop2plot==> 'connect','extRad',...)
% from a specific data run dimensionality (2D or 3D). This is more of a
% convenience function that is used within plotPoreProps.

function [data2Plot] = getData2Plot(file2Analyze,dim,prop2Plot)

    [~,nameNoExt,~] = cellfun(@fileparts,{file2Analyze.name},'UniformOutput',0);
    for i = 1:length(file2Analyze)
        
        currentFile = [file2Analyze(i).folder filesep file2Analyze(i).name];
        field = ['cond' nameNoExt{i}];
        tmp = load(currentFile,'allData');
        switch dim
            case '2D'
                data2Plot.(field) = [tmp.allData.pores2D];

            case '3D'
                data2Plot.(field) = [tmp.allData.pores3D];
            otherwise
                error('Unknown analysis dimension, Please provide correct analysis dimension (2D or 3D)');
        end
        data2Plot.(field) = cell2table({data2Plot.(field).(prop2Plot)}','VariableNames',{prop2Plot});
    end
    
end