%The aim of this code is to save binary tiff files based on a binary
%image and a filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boris Louis (https://github.com/BorisLouis)                             %
% Rafael Camacho Dejay (https://github.com/CamachoDejay)                  %
%                                                                         %
% Website : Rafael Camacho Dejay: https://camachodejay.github.io/         %
%           Boris Louis: https://borislouis.github.io/                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%There is no output to this function. A binary tiff is saved in the current
%matlab folder if just a filename is provided. For specific storage one
%just need to provide tifName as : 'FolderName\tifName'

function BinaryTiff(tifName,BW)
    assert(ischar(tifName),'Filename needs to be a char');
    assert(ismember(length(size(BW)),[2 3]), 'The data you are trying to save has an unexpected dimension');
    test = unique(BW);
    assert(length(test)>2,'The input data is not binary');

    t = Tiff(tifName, 'w');
    setTag(t,'ImageLength',size(BW,1))
    setTag(t,'ImageWidth',size(BW,2))
    setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
    setTag(t,'BitsPerSample',1)
    setTag(t,'SamplesPerPixel',1)
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
    t.write(BW(:,:,1))

    for i = 2:size(BW,3)
        t.writeDirectory
        setTag(t,'ImageLength',size(BW,1))
        setTag(t,'ImageWidth',size(BW,2))
        setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
        setTag(t,'BitsPerSample',1)
        setTag(t,'SamplesPerPixel',1)
        setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
        t.write(BW(:,:,i))
    end
    t.close 
end