function [M,W, MFlow, WFlow] = caMotionClustering(OptSavePath,ResultsPath,SeqNum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  M: Motion Vector Field
%  W: The Weight of each dominant direction at each point. The size is the
%     same as the one of M.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileList = dir(OptSavePath);
cntNum = 0;

%start parallel computing

for i = 1:1:length(fileList)
%for i = 1:1:250
    flowName = fileList(i).name;
    
    if (~strcmp(flowName, '.') && ~strcmp(flowName, '..'))
       [~,name] = fileparts(flowName);
       flowFilePath = fullfile(OptSavePath,[name,'.flo']);
    
       if (exist(flowFilePath, 'file'))
           uv = readFlowFile(flowFilePath);
           
           cntNum = cntNum + 1;
           
           [sizeM, sizeN, ~] = size(uv);
           if sizeN > 1600
              uvS(:,:,1) = imresize(uv(:,:,1),0.25,'nearest');
              uvS(:,:,2) = imresize(uv(:,:,2),0.25,'nearest');
           elseif sizeN > 800
              uvS(:,:,1) = imresize(uv(:,:,1),0.5,'nearest');
              uvS(:,:,2) = imresize(uv(:,:,2),0.5,'nearest');
           else
              uvS(:,:,1) = uv(:,:,1);
              uvS(:,:,2) = uv(:,:,2); 
           end

           %discard small magnitude optical flow
           uvS = caFlowLargeFilter(uvS,0.2);
 
           %normalize optical flow
           uvNorm = caFlowNormlize(uvS);
           
           xFlow{cntNum} = uvNorm(:,:,1);
           yFlow{cntNum} = uvNorm(:,:,2);

       end
    end
end


TotalNumOfFrame = cntNum;
[sizeM, sizeN, ~] = size(uvS);
M{sizeM,sizeN} = 0;
W{sizeM,sizeN} = 0;
MFlow{1} = zeros(sizeM,sizeN,2);
WFlow{1} = zeros(sizeM,sizeN);

for indexR = 1:1:sizeM
    for indexC = 1:1:sizeN
       
        cntNZ = 0;
        for frameNum = 1:1:cntNum                
            vX = xFlow{frameNum}(indexR,indexC);
            vY = yFlow{frameNum}(indexR,indexC);
            if vX~=0 && vY~=0 %discard zero optical flow
              cntNZ = cntNZ + 1;
              vFlow(1,cntNZ) = vX;
              vFlow(2,cntNZ) = vY;   
            end
        end 
        
        if cntNZ > 0
             bandwidth = 0.3;
             [clustCent,~,clustMembsCell] = MeanShiftFlowCluster(vFlow,bandwidth);
             numClust = length(clustMembsCell);
             for i = 1:1:numClust
                 numOfMembs(i) = length(clustMembsCell{i,1});
             end

             [sNum, indexNum] = sort(numOfMembs, 'descend');
             sClustCent = clustCent(:,indexNum);

             M{indexR, indexC} = sClustCent;
             W{indexR, indexC} = sNum;

             for i = 1:1:numClust
                 
                 if i > length(MFlow)
                    MFlow{i} = zeros(sizeM,sizeN,2);
                    WFlow{i} = zeros(sizeM,sizeN);
                 end
                 
                 MFlow{i}(indexR, indexC,1) = sClustCent(1,i) * sNum(i)/cntNum;
                 MFlow{i}(indexR, indexC,2) = sClustCent(2,i) * sNum(i)/cntNum;

                 WFlow{i}(indexR, indexC) = sNum(i)/cntNum;
             end

    %         MSum(indexR, indexC,1) = (sClustCent(1,:) * sNum')/cntNum;
    %         MSum(indexR, indexC,2) = (sClustCent(2,:) * sNum')/cntNum;
    % 
    %         WSum(indexR, indexC) = sum(sNum)/cntNum;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Show the clustering results
    %         figure(10),clf,hold on
    %         cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
    %         for k = 1:min(numClust,length(cVec))
    %             myMembers = clustMembsCell{k};
    %             myClustCen = clustCent(:,k);
    %             plot(x(1,myMembers),x(2,myMembers),[cVec(k) '.'])
    %             plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
    %         end
    %         title(['no shifting, numClust:' int2str(numClust)])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear clustCent;
            clear numOfMembs;
            clear vFlow;
     
            
        end
        
        
        %clc;
        processPercent = 100*((indexR - 1)*sizeN + indexC)/(sizeM*sizeN);
        str = sprintf('Seq%02d---Motion Clustering --- %0.2f%s',SeqNum,processPercent,'%');
        disp(str);

    end
end


%%%save flow field and the corresponding weights
MFlowName = 'MFlow.mat';
WFlowName = 'WFlow.mat';

MFlowSavePath = fullfile(ResultsPath,MFlowName);
WFlowSavePath = fullfile(ResultsPath,WFlowName);

save(MFlowSavePath, 'MFlow');
save(WFlowSavePath, 'WFlow');
%%%



