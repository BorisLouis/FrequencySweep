function trackRes = convertData2TrackRes(allData,nParticles)
  %need to improve this function
  %store only the path in main like in trackRes
  %then have all traces in one together with their respective filenames
  trackRes = struct();
  trackRes.path = allData(1).path;
  traces = cell(nParticles*size(allData,2),2);
    for i =1: size(allData,2)
        nFrames = length(allData(i).traces);
        for j = 1: nParticles
            tabData = table(zeros(nFrames,1),zeros(nFrames,1),zeros(nFrames,1),...
                zeros(nFrames,1),'VariableNames',{'row','col','z','t'});

            tabData.col = allData(i).traces(:,1,j);
            tabData.row = allData(i).traces(:,2,j);
            if size(allData(i).traces,2) <3
                tabData.z = nan(size(allData(i).traces(:,1,j)));
                
            else
                tabData.z = allData(i).traces(:,3,j);
            end
            
            tabData.t = (1:nFrames)';

            traces{(i-1)*nParticles+j,1} = tabData;
            traces{(i-1)*nParticles+j,2} = i;
        end
        
        
    end
    trackRes.traces = traces;


end