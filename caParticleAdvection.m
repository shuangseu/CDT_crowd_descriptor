
function[subMotionIndexTotal] = caParticleAdvection(SeqSourcePath,ResultsPath,MFlow,TimeLength,para)

    
cntNum = 0;

 [m, n] = size(MFlow{1}(:,:,1));
 flowNum = length(MFlow);
 iterNum = flowNum;
 
 correlate = zeros(m,n,flowNum);

%Paticle Advection for calculating flow map
for i = 1:1:TimeLength
     
           cntNum = cntNum + 1;
           if cntNum == 1 
                
                [Xmesh Ymesh] = meshgrid(1:n, 1:m);
   
                for iter = 1:1:flowNum
                    xFlowMap{iter,cntNum} = single(Xmesh);
                    yFlowMap{iter,cntNum} = single(Ymesh); 
                end
                
                  uvPre = MFlow;
                  uvNext = uvPre;
           else
              
               %choose the next velocity 
               for iter = 1:1:iterNum  
                   
                   for indexFlow = 1:1:flowNum

                       MFlowCurrent{indexFlow}(:,:,1) = interp2(Xmesh,Ymesh,MFlow{indexFlow}(:,:,1),xFlowMap{iter,cntNum-1},yFlowMap{iter,cntNum-1},'nearest',0);
                       MFlowCurrent{indexFlow}(:,:,2) = interp2(Xmesh,Ymesh,MFlow{indexFlow}(:,:,2),xFlowMap{iter,cntNum-1},yFlowMap{iter,cntNum-1},'nearest',0);

                       correlate(:,:,indexFlow) = uvPre{iter}(:,:,1).* MFlowCurrent{indexFlow}(:,:,1) + uvPre{iter}(:,:,2).* MFlowCurrent{indexFlow}(:,:,2);
                   end

                   [correlateOrder, indexOrder] = sort(correlate(:,:,:), 3, 'descend');
                   indexLargeCorr = correlateOrder(:,:,1) > 0.6;

                   uvNext{iter}(:,:,1) = zeros(m,n);
                   uvNext{iter}(:,:,2) = zeros(m,n);
                    for indexFlow = 1:1:flowNum
                        uvNext{iter}(:,:,1) = uvNext{iter}(:,:,1) + MFlowCurrent{indexFlow}(:,:,1).*(indexOrder(:,:,1) == indexFlow).*indexLargeCorr;
                        uvNext{iter}(:,:,2) = uvNext{iter}(:,:,2) + MFlowCurrent{indexFlow}(:,:,2).*(indexOrder(:,:,1) == indexFlow).*indexLargeCorr;       
                    end

%                     uvNext{iter}(:,:,1) = uvNext{iter}(:,:,1).*indexLargeCorr;
%                     uvNext{iter}(:,:,2) = uvNext{iter}(:,:,2).*indexLargeCorr;
                    
                    indexFlowMap{iter,cntNum} = indexLargeCorr.* indexOrder(:,:,1);
                    
                   if para.reverseAdvect
                    xFlowMap{iter,cntNum} = (xFlowMap{iter,cntNum-1} - uvNext{iter}(:,:,1));
                    yFlowMap{iter,cntNum} = (yFlowMap{iter,cntNum-1} - uvNext{iter}(:,:,2));
                  else
                    xFlowMap{iter,cntNum} = (xFlowMap{iter,cntNum-1} + uvNext{iter}(:,:,1));
                    yFlowMap{iter,cntNum} = (yFlowMap{iter,cntNum-1} + uvNext{iter}(:,:,2)); 
                   end
               end
               
               uvPre = uvNext;
               
               uvSum = zeros(m,n,2);
               for iter = 1:1:iterNum
                   uvSum = uvSum + abs(uvNext{iter});
               end
               magSum = uvSum(:,:,1) + uvSum(:,:,2);
               if isempty(find(magSum >0,1)) %if the uvNext is equal to zero,then stop the advection
                   break;
               end

           end
            
%         processPercent = 100*i/TimeLength;
%         str = sprintf('Paticle Advection --- %0.2f%s',processPercent,'%');
%         disp(str); 
          str = sprintf('Particle Advection --- iterNum:%d',cntNum);
          disp(str);
   end


FrameNumber = cntNum;

%display the trajectory of particle
if para.display
    %load an example image
    imFileList = dir(SeqSourcePath);
    for i = 1:1:length(imFileList)
        imFileName = imFileList(i).name;
    
        if (~strcmp(imFileName, '.') && ~strcmp(imFileName, '..'))
            imFilePath = fullfile(SeqSourcePath, imFileName);
            im = imread(imFilePath);
            imR = imresize(im,[m,n],'nearest');
            if ~isempty(imR)
                break;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot the trajectory of particle
    x = zeros(1,FrameNumber);
    y = zeros(1,FrameNumber);

    step = 20;

    handleFig = figure();
    imshow(imR);
    hold on;

  for iter = 1:1:iterNum
    for i = 1:step:m
      for j = 1:step:n
       
        for k = 1:FrameNumber

            x(k) = (xFlowMap{iter,k}(i,j));
            y(k) = (yFlowMap{iter,k}(i,j));

        end 

        colorID = mod(iter + i*n + j, 7);

        switch colorID
            case 0
                 sColor = 'r';
            case 1
                 sColor = 'g';
            case 2
                 sColor = 'b';
            case 3
                 sColor = 'c';
            case 4
                 sColor = 'm';
            case 5
                 sColor = 'y';
            case 6
                 sColor = 'w';
        end

        plot(x,y,'.-','LineWidth',4,'Color',sColor);
        hold on;

      end
    end 
  end
    
    FigName = 'ParticleAdvection.png';
    FigSavePath = fullfile(ResultsPath, FigName);
    print(handleFig,'-dpng',FigSavePath);
end


%Density Map%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 densityMap = zeros(m,n);
for iter = 1:1:iterNum
    xParticleMap = round(xFlowMap{iter,FrameNumber});
    yParticleMap = round(yFlowMap{iter,FrameNumber});
    xParticleMap(xParticleMap > n) = n; xParticleMap(xParticleMap < 1) = 1; 
    yParticleMap(yParticleMap > m) = m; yParticleMap(yParticleMap < 1) = 1;
    
    for i = 1:m
        for j = 1:n
            densityMap(yParticleMap(i,j),xParticleMap(i,j)) = densityMap(yParticleMap(i,j),xParticleMap(i,j)) + 1;        
        end
    end
end

%Gaussian filter
hG = fspecial('gaussian',[11 11],1);
densityMapGaussian = imfilter(densityMap,hG);

%show densityMap and GaussianDensityMap
if para.display
    figure();
    imagesc(densityMap);

    handleFig = figure();
    imagesc(densityMapGaussian);

    FigName = 'densityMapGaussian.png';
    FigSavePath = fullfile(ResultsPath, FigName);
    print(handleFig,'-dpng',FigSavePath);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Local Peaks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MeanShift Cluster
MaxDensity = max(max(densityMapGaussian));
thresholdForLocalPeaks = MaxDensity * para.threForLocalPeaks;

%thresholdForLocalPeaks =  para.threForLocalPeaks;

bandWidth = 4;

[clustCent,point2cluster,clustMembsCell] = MeanShiftDensityCluster(densityMapGaussian,thresholdForLocalPeaks, bandWidth);

if para.display
    handleFig = figure();
    imagesc(densityMapGaussian);
    hold on;

    LocalPeaks = round(clustCent);
    numClust = length(clustMembsCell);
    for k = 1:numClust
        plot(LocalPeaks(2,k),LocalPeaks(1,k),'--rs','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','g', 'MarkerSize',10)
    end
    title(['Local Peaks, Number:' int2str(numClust)]);
end

%Accumulation points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[accPts,accPoint2cluster,accClustMembsCell] = AccumulationPts(densityMapGaussian, LocalPeaks,  para.widthForAccPts);

weight4Sig = zeros(1,length(accClustMembsCell));
for i = 1:length(accClustMembsCell)
    index =LocalPeaks(:,accClustMembsCell{i,1});
    
    maxRow = max(index(1,:));
    minRow = min(index(1,:));
    maxCol = max(index(2,:));
    minCol = min(index(2,:));

    RectMaxRow = min(m, maxRow + 10);
    RectMinRow = max(1, minRow - 10);
    RectMaxCol = min(n, maxCol + 10);
    RectMinCol = max(1, minCol - 10);
    
    MaskR = zeros(m,n);
    for ii = RectMinRow:RectMaxRow
        for jj= RectMinCol:RectMaxCol
            MaskR(ii,jj) = 1;
        end
    end
    
    densityROI = densityMapGaussian.*MaskR;
    weight4Sig(i) = sum(sum(densityROI));
    
%     [indexR,indexC] = size(index);
%     for j = 1:indexC
%         weight4Sig(i) = weight4Sig(i) + densityMapGaussian(index(1,j),index(2,j));
%     end
end

thresholdForAccPts =  para.threForAccPts;
indexAccPts = find(weight4Sig > thresholdForAccPts);
selectAccPts = zeros(2,length(indexAccPts));

if para.display
    
    for i = 1:length(indexAccPts)
    selectAccPts(:,i) = round(accPts(:,indexAccPts(i)));
    plot(selectAccPts(2,i), selectAccPts(1,i),'--rs','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',10);  
    end
    
    FigName = 'LocalPeaks.png';
    FigSavePath = fullfile(ResultsPath, FigName);
    print(handleFig,'-dpng',FigSavePath);

end

%tajectoryDirection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumSelectAccPts = length(indexAccPts);
for i = 1:1:NumSelectAccPts %compute trajectory for each accumulation point.
    
    clear trejectoryHis;
    for iter = 1:1:iterNum  %comptue trajectory for iterNum motion layers
    
    LocalPeaksOfAccPts = LocalPeaks(:,accClustMembsCell{indexAccPts(i),1});
    
    [theta, MaskR,Mask4NearAccPts, Mask4LongTrajectory, trajectoryLen] = caTrajectoryDirectionHist(xFlowMap(iter,:),yFlowMap(iter,:),densityMapGaussian, LocalPeaksOfAccPts,accPts,10, para.threForLongTrajectory);
    
%   handleFig = figure();
    bins = 0:10:360;
    trajectoryHist{iter} = hist(theta(Mask4LongTrajectory),bins);
   
%     FigName = sprintf('trajectoryHist%d.png',i);
%     FigSavePath = fullfile(ResultsPath, FigName);
%     print(handleFig,'-dpng',FigSavePath);
    
%     handleFig = figure();
%     imagesc(MaskR);
%             
%     FigName = sprintf('MaskRegionNearAccPts%d.png',i);
%     FigSavePath = fullfile(ResultsPath, FigName);
%     print(handleFig,'-dpng',FigSavePath);
%     
%     if ~isempty(find(Mask4LongTrajectory>0))
%     direcData = theta(Mask4LongTrajectory>0);
%     direcDataPts = zeros(2,length(direcData));
%     direcDataPts(1,:) = direcData;
%     direcDataPts(2,:) = 0;
%     direcBW = 20;
%     [direcClustCent,direcData2cluster,direcCluster2dataCell] = MeanShiftDirectionCluster(direcDataPts,direcBW);
%     end
    
    %Mask for Path
    %find submotion flow
    Mask4Path = zeros(m,n);
    [xIndex, yIndex] = find(Mask4LongTrajectory>0);

     x = zeros(1,FrameNumber);
     y = zeros(1,FrameNumber);

     subMotionIndexTotal{i,iter} = zeros(m,n,flowNum);

    step = 20;

%     if para.display
%         figure();
%         imshow(im);
%         hold on;
%     end
            
%     for ii = 1:step:length(xIndex)
% 
%         for k = 1:FrameNumber
% 
%             x(k) = round(xFlowMap{k}(xIndex(ii),yIndex(ii)));
%             y(k) = round(yFlowMap{k}(xIndex(ii),yIndex(ii)));
%             Mask4Path(y(k),x(k)) = 1;
%         end 
% 
%         if para.display
%             colorID = mod(ii, 7);
%             switch colorID
%                 case 0
%                      sColor = 'r';
%                 case 1
%                      sColor = 'g';
%                 case 2
%                      sColor = 'b';
%                 case 3
%                      sColor = 'c';
%                 case 4
%                      sColor = 'm';
%                 case 5
%                      sColor = 'y';
%                 case 6
%                      sColor = 'w';
%             end
% 
%             plot(x,y,'.-','LineWidth',4,'Color',sColor);
%             hold on;
%         end
% 
%     end

    for ii = 1:2:length(xIndex)
        
        x(1) = round(xFlowMap{iter,1}(xIndex(ii),yIndex(ii)));
        y(1) = round(yFlowMap{iter,1}(xIndex(ii),yIndex(ii)));
        
    for k = 2:FrameNumber

        x(k) = round(xFlowMap{iter,k}(xIndex(ii),yIndex(ii)));
        y(k) = round(yFlowMap{iter,k}(xIndex(ii),yIndex(ii)));
        
        if x(k)>n
            x(k) = n;
        end
        
        if x(k)<1
            x(k) = 1;
        end
        
        if y(k)>m
            y(k) = m;
        end
        
        if y(k)<1
            y(k) = 1;
        end
        
        Mask4Path(y(k),x(k)) = 1;
        
        %compute the submotion index 
        indexFlow = indexFlowMap{iter,k}(xIndex(ii),yIndex(ii));
        if indexFlow>0
           subMotionIndexTotal{i,iter}(y(k-1),x(k-1),indexFlow) = subMotionIndexTotal{i,iter}(y(k-1),x(k-1),indexFlow) + 1;
        end  
        
        
    end 
    end
        
    se1=strel('disk',5);
    Mask4PathROI{i} = imdilate(Mask4Path,se1);
%     if para.display
%         handleFig = figure();
%         imagesc(Mask4PathROI{i});
%         
%         FigName =sprintf('Mask4PathROI%d.png',i);
%         FigSavePath = fullfile(ResultsPath, FigName);
%         print(handleFig,'-dpng',FigSavePath);
%     end
    end 
    
    clear histSum;
    for iter = 1:1:iterNum
        if iter == 1
           histSum = trajectoryHist{iter};
        else
           histSum = histSum + trajectoryHist{iter};
        end  
    end
    
    bins = 0:10:360;
    handleFig = figure();
    bar(bins, histSum);
    FigName = sprintf('trajectoryHist%d.png',i);
    FigSavePath = fullfile(ResultsPath, FigName);
    print(handleFig,'-dpng',FigSavePath);
  
end

accPtsName = 'accPts.mat';
mask4PathROIName = 'mask4PathROI.mat';

if  para.reverseAdvect
    paraName = 'R_para.mat';
else
    paraName = 'para.mat';
end

accPtsSavePath = fullfile(ResultsPath,accPtsName);
maskSavePath = fullfile(ResultsPath,mask4PathROIName);
paraSavePath = fullfile(ResultsPath,paraName);
subMotionIndexSavePath = fullfile(ResultsPath,'subMotionIndexTotal.mat');

save(accPtsSavePath, 'selectAccPts');
save(maskSavePath, 'Mask4PathROI'); 
save(paraSavePath, 'para');
save(subMotionIndexSavePath,'subMotionIndexTotal');

str = sprintf('Paticle Advection --- End',cntNum);
disp(str);







    




