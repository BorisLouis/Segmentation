

%% Load hot

p2file      = 'D:\Documents\Unif\PhD\2022-Data\07 - July\Fiber analysis\Hot\Fr1.tif';
warning('off','all')
fileInfo    = Load.Movie.tif.getinfo(p2file);

warning('on','all')
tNframes = fileInfo.Frame_n;

%loop through the frames of the current stack
nIM = tNframes;
% waitbar(j/nImStacks,h,hMessage);
frames = 1:fileInfo.Frame_n;
hotData = Load.Movie.tif.getframes(p2file,frames);

%% Load Cold
p2file      = 'D:\Documents\Unif\PhD\2022-Data\07 - July\Fiber analysis\Cold\Fr1.tif';
warning('off','all')
fileInfo    = Load.Movie.tif.getinfo(p2file);

warning('on','all')
tNframes = fileInfo.Frame_n;

%loop through the frames of the current stack
nIM = tNframes;
% waitbar(j/nImStacks,h,hMessage);
frames = 1:fileInfo.Frame_n;
coldData = Load.Movie.tif.getframes(p2file,frames);

%% Correlation analysis
[corrMat] = normxcorr2(hotData(:,:,1),hotData(:,:,1));
[a1,b1] = Calc.radial_profile(corrMat,1);

[corrMat] = normxcorr2(coldData(:,:,1),coldData(:,:,1));
[a2,b2] = Calc.radial_profile(corrMat,1);

figure
plot(a2+1,b2/sum(b2))
hold on
plot(a1+1,b1/sum(b1))
set(gca,'XScale','log')
legend({'Cold','Hot'})

%% without preProcess

 %test fft
[hotFFT] = fftshift(fft2(double(hotData)));

[hotA,hotB] = Calc.radial_profile(abs(hotFFT),1);



[coldFFT] = fftshift(fft2(double(coldData)));
[coldA,coldB] = Calc.radial_profile(abs(coldFFT),1);


figure
hold on
plot(coldA,(coldB-19237)./(sum((coldB-19237))))
plot(hotA,(hotB-26000)./(sum((hotB-26000))))
% plot(1./coldA,coldB)
% plot(1./hotA,hotB)
hold on

set(gca,'XScale','log')
set(gca,'YScale','log')
axis square
box on


%% with prePorcess
tmpData =double(imgaussfilt(hotData));

figure
test2 = fibermetric(tmpData,'StructureSensitivity',10);
subplot(1,2,1)
imagesc(test2)
axis image
colormap('hot')

tmpColdData =double(imgaussfilt(coldData));
test3 = fibermetric(tmpColdData,'StructureSensitivity',10);
subplot(1,2,2)
imagesc(test3)
axis image
colormap('hot')

[BWCold] = imbinarize(test3);
[BWHot]  = imbinarize(test2);
figure
subplot(1,2,1)
imagesc(BWCold)

subplot(1,2,2)
imagesc(BWHot)



[~,~,bwCold] = imSegmentation.segmentStack(test3,'threshold',0.6,'connectivity',501);
[~,~,bwHot]  = imSegmentation.segmentStack(test2,'threshold',0.6,'connectivity',501);

% clean bw
bwCold = bwareaopen(bwCold,10);
bwHot = bwareaopen(bwHot,10);

SE = strel('disk',5);
bwCold = imclose(bwCold,SE);
bwHot = imclose(bwHot,SE);



% bwCold = imfill(bwCold,'holes');
% bwHot  = imfill(bwHot,'holes');

figure
subplot(1,2,1)
imagesc(bwCold)

subplot(1,2,2)
imagesc(bwHot)


%% gradient
figure

[Gx,Gy] = imgradientxy(tmpData);
[Gmag,Gdir] = imgradient(Gx,Gy);
subplot(2,2,1)
imagesc(Gmag)
subplot(2,2,2)
imagesc(Gdir)

[Gx,Gy] = imgradientxy(tmpColdData);
[Gmag,Gdir] = imgradient(Gx,Gy);
subplot(2,2,3)
imagesc(Gmag)
subplot(2,2,4)
imagesc(Gdir)



%%
[hotFiltFFT] = fftshift(fft2(double(imgaussfilt(hotData,3))));
[hotFiltA,hotFiltB] = Calc.radial_profile(abs(hotFiltFFT),1);

[coldFiltFFT] = fftshift(fft2(double(imgaussfilt(coldData,3))));
[coldFiltA,coldFiltB] = Calc.radial_profile(abs(coldFiltFFT),1);

figure
hold on
plot(1./hotFiltA,hotFiltB)
hold on
plot(1./coldFiltA,coldFiltB)
set(gca,'XScale','log')
set(gca,'YScale','log')
axis square
box on
