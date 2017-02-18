function [subMergeMotionField] = caExtractSubMotion(subMotionIndexTotal,MFlowDeno,ResultsPath)

str = sprintf('Extract sub-Motion Field --- Start.');
disp(str);

[sizeM,sizeN] = size(subMotionIndexTotal);
for cntSubMIT = 1:1:sizeM
   subMotionIndex = subMotionIndexTotal(cntSubMIT,:);
    
%     %compute the submotion
    numOfMotion = length(subMotionIndex);
    for indexMotion = 1:1:numOfMotion
  
    [m,n,flowNum] = size(subMotionIndex{indexMotion});
    subMotionField = zeros(m,n,2);
    subMotionMask = zeros(m,n);
    
    [subMotionOrder, indexOrder] = sort(subMotionIndex{indexMotion}(:,:,:),3, 'descend');
    subMotionMask = (subMotionOrder(:,:,1) > 0);
    
    for indexFlow = 1:1:flowNum  
        subMotionField(:,:,1) = subMotionField(:,:,1) + MFlowDeno{indexFlow}(:,:,1).*(indexOrder(:,:,1) == indexFlow).*subMotionMask;
        subMotionField(:,:,2) = subMotionField(:,:,2) + MFlowDeno{indexFlow}(:,:,2).*(indexOrder(:,:,1) == indexFlow).*subMotionMask;
    end
    
%     figure;
%     imagesc(subMotionMask);
%     
%      figure;
%      plotflow(subMotionField);
%     
%     figure;
%     magnitude = abs(subMotionField(:,:,1)) + abs(subMotionField(:,:,2));
%     imagesc(magnitude>0);
    
    %subMotoin medfilter
    subMotionFieldMed(:,:,1) = medfilt2(subMotionField(:,:,1),[5,5]);
    subMotionFieldMed(:,:,2) = medfilt2(subMotionField(:,:,2),[5,5]);
    
    subMotionFieldMed(:,:,1) = medfilt2(subMotionFieldMed(:,:,1),[5,5]);
    subMotionFieldMed(:,:,2) = medfilt2(subMotionFieldMed(:,:,2),[5,5]);
    
    subMotionFieldMed(:,:,1) = medfilt2(subMotionFieldMed(:,:,1),[5,5]);
    subMotionFieldMed(:,:,2) = medfilt2(subMotionFieldMed(:,:,2),[5,5]);
    
    subMotionFieldMed(:,:,1) = medfilt2(subMotionFieldMed(:,:,1),[5,5]);
    subMotionFieldMed(:,:,2) = medfilt2(subMotionFieldMed(:,:,2),[5,5]);
    
%     figure;
%     plotflow(subMotionFieldMed);
% %     
%     figure;
%     magnitude = abs(subMotionFieldMed(:,:,1)) + abs(subMotionFieldMed(:,:,2));
%     imagesc(magnitude>0);
    
    %subMotionFieldMed = caFlowLargeFilter(subMotionFieldMed,0.1);
    subLayer{indexMotion} = caFlowNormlize(subMotionFieldMed);
    
    end
    
    % Merge subMotionFieldLayer
    if length(subLayer)>1
        magnitude = abs(subLayer{1}(:,:,1) + subLayer{1}(:,:,2));
        MaskLayer{1} = (magnitude > 0);
        
        clear subMerge;
        subMerge{1} = subLayer{1};
        for indexLayer = 2:1:length(subLayer)
            magnitude = abs(subMerge{1}(:,:,1) + subMerge{1}(:,:,2));
            MaskLayer{1} = (magnitude > 0);
            
            magnitude = abs(subLayer{indexLayer}(:,:,1) + subLayer{indexLayer}(:,:,2));
            MaskLayer{2} = (magnitude > 0);
            
            sumMask = MaskLayer{1} + MaskLayer{2};
            if length(find(sumMask == 2)) % there are intersection points 
               intersecMask = (sumMask == 2);
               correlate = subMerge{1}(:,:,1).*subLayer{indexLayer}(:,:,1) + subMerge{1}(:,:,2).*subLayer{indexLayer}(:,:,2);
               correForIntersec = correlate.*intersecMask;
               lowCorreMask = (correForIntersec < 0.95).*intersecMask;
               if length(find(lowCorreMask==1)) > 400 %check the correlation of intersection points
                   cntTemp = length(subMerge) + 1;
                   subMerge{cntTemp} = subLayer{indexLayer};
               else
                   %only add the no intersection points
                   MaskForMerge(:,:,1) = 1 - intersecMask;
                   MaskForMerge(:,:,2) = 1 - intersecMask;
                   subMerge{1} = caFlowNormlize(subMerge{1} + subLayer{indexLayer}.* MaskForMerge);
               end
            else %there is no intersection
                subMerge{1} = caFlowNormlize(subMerge{1} + subLayer{indexLayer});
            end
 
        end
   
    else
        subMerge = subLayer;
    end
        

    handleFig = figure;
    for i = 1:1:length(subMerge)    
       subMerge{i}(:,:,1) = medfilt2(subMerge{i}(:,:,1),[10,10]);
       subMerge{i}(:,:,2) = medfilt2(subMerge{i}(:,:,2),[10,10]);
%      subMerge{i}(:,:,1) = medfilt2(subMerge{i}(:,:,1),[10,10]);
%      subMerge{i}(:,:,2) = medfilt2(subMerge{i}(:,:,2),[10,10]);
%       figure;
       plotflow(subMerge{i});
       hold on;   
    end
    hold off;
    
    FigName = sprintf('SubMotionFig-%02d.png',cntSubMIT);
    FigSavePath = fullfile(ResultsPath, FigName);
    print(handleFig,'-dpng',FigSavePath);
    
  subMergeMotionField{cntSubMIT} = subMerge;  
  
end
    
subMotionFieldSavePath = fullfile(ResultsPath, 'subMergeMotionField.mat');
save(subMotionFieldSavePath, 'subMergeMotionField');

str = sprintf('Extract sub-Motion Field --- End.');
disp(str);