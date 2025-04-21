clear 
clc
close all
%% User input
%file info
file.path  = 'E:\Flip\03 - 12 - Freq. Sweep\66nm Andrey Particles - 2';
file.ext   = '.ome.tif';

particleType = 'PS'; %'PS' %'Au'
freqFact = 1000;%factor to multiply the frequency with (
startFreq = 20000;
endFreq = 10;
nFreq = 50;
pxSize = 95;%in nm
info.type = 'normal';%Transmission or normal 
info.checkSync = false; %true if want to check for camera synchronization
info.useSameROI = true;
info.runMethod = 'load';% 'run'
outputFolder = 'Results'; %name of the folder to output the results

%% Loading
tmp = dir(file.path);%make it a directory in matlab (list of all thing contained inside)
tmp = tmp(cell2mat({tmp.isdir}));%only keep the subfolders
tmp = tmp(3:end);%remove the access to parent folders

folder2Mov = [];
for i = 1:size(tmp,1)
    path2File = [tmp(i).folder filesep tmp(i).name];%get path to subfolders
    file2Analyze = Core.Movie.getFileInPath(path2File,'.ome.tif');%get .ome.tif in subfolder
    if size(file2Analyze,1)>1
        file2Analyze = file2Analyze(1);
    end
        folder2Mov = [folder2Mov file2Analyze];%store the info to the file if it found one (otherwise file2Analyze is empty)
end

%check that folder2Mov is not empty
assert(~isempty(folder2Mov),'Error no file was found, check that you put the correct analysis type');

%% make a binary cross
side = 350;

cross = logical(eye(side,side)+imrotate(eye(side,side),90));
crossDetect = imdilate(cross,strel('disk',10));

id = round(size(crossDetect,1)/2)-10:round(size(crossDetect,1)/2)+10;
crossDetect(id,id)=0;
cross = imdilate(cross,strel('disk',4));

% figure 
% imagesc(cross)
freqRange = logspace(log10(startFreq*freqFact),log10(endFreq*freqFact),nFreq);

%% detection of the center of the beads
%preallocate memory for storing data of the different files
allData = struct('fileName',[],'locPos',[]);
allData(size(folder2Mov,2)).locPos = [];

for i =1: size(folder2Mov,2)
    
    f.path = folder2Mov(i).folder;
    f.ext  = '.ome.tif';
    
    currMov = Core.Movie(f,info);%Create Movie Object
    fullStack = currMov.getFrame(1);%extract first frame
    frame = fullStack.Cam1;
    %check if cropping is necessary
    if size(frame,2) > 400
            currMov.cropIm;
            prevROI = currMov.info.ROI;
    else
        %prevROI = [1 ,1, currMov.raw.movInfo.Width,currMov.raw.movInfo.Length];
    end
    %load full stack
    fullStack = currMov.getFrame;
    currentPath = currMov.raw.movInfo.Path;
    %revert the intensity scale
    if strcmpi(info.type,'transmission')
        fullStackIn = imcomplement(fullStack.Cam1);
    else
        fullStackIn = double(fullStack.Cam1);
    end
  %%  detect cross
    meanIm = max(fullStackIn,[],3);
   
    % t = imregcorr(crossDetect,meanIm,'similarity');
    % if abs(t.RotationAngle)>5
    %     t.RotationAngle = -180-t.RotationAngle;
    % end
    % fCross = zeros(size(cross,1)+round(t.Translation(2)),size(cross,2)+round(t.Translation(1)));
    % fCross(round(t.Translation(2)):end-1,round(t.Translation(1)):end-1) = cross;
    % fCross = imrotate(fCross,-t.RotationAngle);
    % fCross(end:size(meanIm,1),end:size(meanIm,2))= 0;

    if strcmpi(particleType,'PS')
        c = normxcorr2(crossDetect,meanIm);
    elseif strcmpi(particleType,'Au')
        c = normxcorr2(~crossDetect,meanIm);
    else
        error('Unkown particle types, PS and Au are the two accepted')
    end
    %delete values too close to the center
    mask = zeros(size(c));
    mask(round(size(c,1)/2)-100:round(size(c,1)/2)+100,round(size(c,2)/2)-100:round(size(c,2)/2)+100) =1;

    c = c.*mask;  

    [ypeak,xpeak] = find(c==max(c(:)));

    yoffSet = round(ypeak-size(crossDetect,1));%cross
    xoffSet = round(xpeak-size(crossDetect,2));

    fCross = zeros(size(cross,1)+yoffSet,size(cross,2)+xoffSet);
    fCross(yoffSet:end-1,xoffSet:end-1) = cross;

    fCross(end:size(meanIm,1),end:size(meanIm,2))= 0;

    figure(20+i)
    imagesc(imfuse(meanIm,fCross))
    title('Alignment of cross with Data')

    %% get Negative DEP cross
    pCross = fCross;
    nCross = imdilate(pCross,strel('disk',20))-pCross;
    
    % figure
     % imagesc(imfuse(pCross,nCross))

    %to normalize by number of pixels
    pNorm = sum(pCross(:));
    nNorm = sum(nCross(:));
    %% detect darkFrames

    fInt = squeeze(sum(sum(fullStackIn,1),2));

    startFrame = find(fInt>min(fInt)*1.1,1,'first');
    
    meanInt = (min(fInt)+max(fInt))/2;
    darkFrame = (fInt>meanInt);
    % figure
    % plot(darkFrame)
    % 
    changepoint =diff(darkFrame);
   % plot(changepoint)
    
    idx2End = find(changepoint==-1)-2;
    idx2start =find(changepoint==1)+2;

    %crop idx
    idx2End=idx2End(1:nFreq);
    idx2start = idx2start(1:nFreq);
    
    int = zeros(length(idx2start),1);
    figure(100)
    for j = 1:length(idx2start)
        % imagesc(imfuse(fullStackIn(:,:,idx2End(j)),nCross))
        % 
        % colorbar
       
        pStart(j) = sum(fullStackIn(:,:,idx2start(j)).*pCross,'all')./pNorm;
        nStart(j) = sum(fullStackIn(:,:,idx2start(j)).*nCross,'all')./nNorm;
        pEnd(j) = sum(fullStackIn(:,:,idx2End(j)).*pCross,'all')./pNorm;
        nEnd(j) = sum(fullStackIn(:,:,idx2End(j)).*nCross,'all')./nNorm;
        % clf
    end

    int= (pEnd-pStart)./(nEnd-nStart);
    pInt   = (pEnd-pStart);
    tInt = pEnd./nEnd;

    % figure
    % plot(pStart)
    % hold on
    % plot(nStart)
    % plot(pEnd)
    % plot(nEnd)

    % figure
    % subplot(1,2,1)
    % scatter(freqRange,int)
    % set(gca,'XScale','log')
    % axis square
    % box on
    % 
    figure
    % subplot(1,2,2)
    scatter(freqRange,pEnd./nEnd,10,'filled')
    set(gca,'XScale','log')
    xlabel('Frequency (Hz)')
    ylabel('Integrated intensity')
    axis square
    box on

    %%




    %store data in allData
    allData(i).pEnd = pEnd;
    allData(i).nEnd = nEnd;
    allData(i).freq = freqRange;
    allData(i).fileName = currentPath;
    allData(i).path = f.path;
    %clear waitbar
   
end

%% save data
meanCurve = zeros(size(pEnd));
for i = 1:length(allData)

    meanCurve = meanCurve + (allData(i).pEnd./allData(i).nEnd);

end

meanCurve = meanCurve/length(allData);

allData(1).meanCurve = meanCurve;

figure
scatter(freqRange,meanCurve,10,'filled')
set(gca,'XScale','log')
xlabel('Frequency (Hz)')
ylabel('Integrated intensity')
axis square
box on

[~,name,~]= fileparts(file.path);

filename = [file.path filesep 'fSweep-' name ];
save(filename,'allData');
h = msgbox('Data succesfully saved');