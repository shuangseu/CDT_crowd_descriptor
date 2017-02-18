function [curlMap,divMap,curlDescriptor, divDescriptorAll,divDescriptorShort, curlInteMap, divInteMap, maskCurlMap, maskDivMap] = caMotionDescription(motionVF)
    curlMap = {};
    divMap ={};
    curlDescriptor = {};
    divDescriptorAll = {};
    divDescriptorShort = {};
    curlInteMap = {};
    divInteMap = {};
    maskCurlMap = {};
    maskDivMap = {};
%%%%%%
% motionVF - motion vector field

%%%%%%
disp('I am computing motion descriptor! Please wait...');

[m,n,z] = size(motionVF);

uv = caFlowNormlize(motionVF);
%normal vector field
NormalUV = zeros(size(uv));
NormalUV(:,:,1) = -uv(:,:,2);
NormalUV(:,:,2) = uv(:,:,1);

%compute the trajectories in tangential and normal directions
disp('-1-.computing the trajectories...');
[trajectoryTan] = caPtsTrajectory(uv,10);
[trajectoryNor] = caPtsTrajectory(NormalUV,2);

numTan = length(trajectoryTan);
numNor = length(trajectoryNor);

if numTan > 0 && numNor > 0

    %compute the intersect points between tangential and normal trajectories
    disp('-2-.computing the intersect points...');
    [intersectPts,intersectPtsTanIndex,intersectPtsNorIndex] = caComputIntersetPts(trajectoryTan,trajectoryNor);

    % figure;
    % caPlotTrajectory(trajectoryNor);
    % plot(intersectPts(:,1),intersectPts(:,2),'o','color','r');

    %name the index by using consecutive numbers. The same index share the same
    %number.
    disp('-3-.sorting the trajectories...');
    [intersectPtsTanIndex2] = caIndexOrderNum(intersectPtsTanIndex);
    [intersectPtsNorIndex2] = caIndexOrderNum(intersectPtsNorIndex.').';

    %sort the intersect points based on the order number.
    [sortColNonZero,sortRowNonZero] = caSortRowColInterPts(intersectPtsTanIndex2,intersectPtsNorIndex2);

    [sortColNonZeroCut] = caSortColIndexCut(intersectPtsTanIndex2,sortColNonZero);
    % caShowSortIndex(intersectPtsTanIndex2,sortColNonZero);
    % caShowSortIndex(intersectPtsNorIndex2.',sortRowNonZero);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('-4-.computing the CD descriptor...');
    clear curlMap;
    clear divMap;
    clear curlDescripotr;
    clear divDescriptorAll;
    clear divDescriptorShort;
    clear curlInteMap;
    clear divInteMap;
    curlMap = {};
    divMap ={};
    curlDescriptor = {};
    divDescriptorAll = {};
    divDescriptorShort = {};
    curlInteMap = {};
    divInteMap = {};
    maskCurlMap = {};
    maskDivMap = {};
    for cnt = 1:1:length(sortColNonZero)
    [curlMap{cnt},divMap{cnt},curlDescriptor{cnt}, divDescriptorAll{cnt},divDescriptorShort{cnt},curlInteMap{cnt}, divInteMap{cnt},maskCurlMap{cnt},maskDivMap{cnt}] = caCDIntegrationDescriptor(uv,sortRowNonZero{cnt},sortColNonZero{cnt},sortColNonZeroCut{cnt},trajectoryTan,trajectoryNor,intersectPts);
    end

else
    disp('Trajectory is too short.');
end

disp('--- END ---');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%functions for motion description
function [maskForPts] = caMotionBeginEndPts(uv,flag)
[m,n,z] = size(uv);
magnitudeMatrix = abs(uv(:,:,1)) + abs(uv(:,:,2));
maskForNonZero = magnitudeMatrix > 0;

[Xmesh Ymesh] = meshgrid(1:n, 1:m);
xFlowMap{1} = single(Xmesh);
yFlowMap{1} = single(Ymesh); 

%flag = 1, forward to find end pts
%flag = -1, backward to find begin pts
xFlowMap{2} = xFlowMap{1} + flag * interp2(Xmesh,Ymesh,uv(:,:,1),xFlowMap{1},yFlowMap{1},'nearest',0);
yFlowMap{2} = yFlowMap{1} + flag * interp2(Xmesh,Ymesh,uv(:,:,2),xFlowMap{1},yFlowMap{1},'nearest',0);

uvNext(:,:,1) =  interp2(Xmesh,Ymesh,uv(:,:,1),xFlowMap{2},yFlowMap{2},'nearest',0);
uvNext(:,:,2) =  interp2(Xmesh,Ymesh,uv(:,:,2),xFlowMap{2},yFlowMap{2},'nearest',0);

magnitudeUvNext = abs(uvNext(:,:,1)) + abs(uvNext(:,:,2));
maskForPts = (magnitudeUvNext == 0).* maskForNonZero;
        
end

function [trajectoryMatrix] = caPtsTrajectory(uv,Th4Length)
trajectoryMatrix = {};

[m,n,z] = size(uv);

[maskForTanBeginPts] = caMotionBeginEndPts(uv,-1);
[maskForTanEndPts] = caMotionBeginEndPts(uv,1);

% for i = 1:4:m
%    maskForTanBeginPts(i,:) = 0; 
% end
% for j = 1:4:n
%    maskForTanEndPts(:,j) = 0;
% end

[irB,icB] = find(maskForTanBeginPts == 1);
[irE,icE] = find(maskForTanEndPts == 1);

endPts = [irE,icE];
maskVisitedEndPts = zeros(1,length(irE));
% clear irB;
% clear icB;
% 
% irB = irE;
% icB = icE;

lengthBPts = length(irB);

[Xmesh Ymesh] = meshgrid(1:n, 1:m);
traNum = 1;
for cnt = 1:1:lengthBPts
    iter  = 1;
    while 1
        if iter == 1
           traForward{cnt}(iter,1) = icB(cnt);
           traForward{cnt}(iter,2) = irB(cnt); 
        else

           uvNext(:,:,1) =  interp2(Xmesh,Ymesh,uv(:,:,1),traForward{cnt}(iter-1,1), traForward{cnt}(iter-1,2),'nearest',0);
           uvNext(:,:,2) =  interp2(Xmesh,Ymesh,uv(:,:,2),traForward{cnt}(iter-1,1), traForward{cnt}(iter-1,2),'nearest',0);
           
           if iter == 2
               uvOld = uvNext;
           end
           
           correlation = uvOld(:,:,1)*uvNext(:,:,1) + uvOld(:,:,2)*uvNext(:,:,2);

%            if uvNext(:,:,1) == 0 && uvNext(:,:,2) == 0
           if correlation < 0.17 % 80 degrees
               break;
           else
              traForward{cnt}(iter,1) = (traForward{cnt}(iter -1,1) + uvNext(:,:,1));
              traForward{cnt}(iter,2) = (traForward{cnt}(iter -1,2) + uvNext(:,:,2)); 
              
              if iter>2 %局部来回循环
                  if(traForward{cnt}(iter,:) == traForward{cnt}(iter -2,:))
                   break;
                  end
              end
             
           end
           uvOld = uvNext;
        end
     iter = iter + 1;
     
     if iter > (2*(m+n))
         break;
     end
%      str = sprintf('cnt =  %03d, iter = %04d',cnt,iter);
%      disp(str);
    end
    
    if iter > Th4Length %ignore the short tracklets
       trajectoryMatrix{traNum} = round(traForward{cnt});
      
       %%%%find the unvisted end points
       [~,~,iE] = intersect(trajectoryMatrix{traNum},endPts,'rows');
       if ~isempty(iE)
           maskVisitedEndPts(1,iE) = 1;
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       traNum = traNum + 1;
    end  
    
end

%%%%%backward particle advection
unVisitedEndPts = endPts(maskVisitedEndPts == 0,:);
if ~isempty(unVisitedEndPts)
    
    lengthEPts = length(unVisitedEndPts(:,1));

    [Xmesh Ymesh] = meshgrid(1:n, 1:m);
%     traNum = 1;
    for cnt = 1:1:lengthEPts
        iter  = 1;
        while 1
            if iter == 1
               traForward{cnt}(iter,1) = unVisitedEndPts(cnt,1);
               traForward{cnt}(iter,2) = unVisitedEndPts(cnt,2); 
               
            else

               uvNext(:,:,1) =  interp2(Xmesh,Ymesh,uv(:,:,1),traForward{cnt}(iter-1,1), traForward{cnt}(iter-1,2),'nearest',0);
               uvNext(:,:,2) =  interp2(Xmesh,Ymesh,uv(:,:,2),traForward{cnt}(iter-1,1), traForward{cnt}(iter-1,2),'nearest',0);
             
               if iter == 2
                  uvOld = uvNext;
               end
           
              correlation = uvOld(:,:,1)*uvNext(:,:,1) + uvOld(:,:,2)*uvNext(:,:,2);

%              if uvNext(:,:,1) == 0 && uvNext(:,:,2) == 0
               if correlation < 0.17 % 80 degrees
                   break;
               else
                  traForward{cnt}(iter,1) = (traForward{cnt}(iter -1,1) - uvNext(:,:,1)); % minus means backward advection
                  traForward{cnt}(iter,2) = (traForward{cnt}(iter -1,2) - uvNext(:,:,2)); 

                  if iter>2 %局部来回循环
                      if(traForward{cnt}(iter,:) == traForward{cnt}(iter -2,:))
                       break;
                      end
                  end

               end
               uvOld = uvNext;
            end
         iter = iter + 1;
        end

        if iter > Th4Length %ignore the short tracklets
           trajectoryMatrix{traNum} = fliplr(round(traForward{cnt}));
           traNum = traNum + 1;
        end  
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if plotFlag
%     figure;
%     for cnt = 1:1:length(traForward)
%         
%          colorID = mod(cnt, 6);
% 
%         switch colorID
%             case 0
%                  sColor = 'r';
%             case 1
%                  sColor = 'g';
%             case 2
%                  sColor = 'b';
%             case 3
%                  sColor = 'c';
%             case 4
%                  sColor = 'm';
%             case 5
%                  sColor = 'y';
%             case 6
%                  sColor = 'w';
%         end
% 
%         plot(round(traForward{cnt}(:,1)),round(traForward{cnt}(:,2)),'.-','LineWidth',4,'Color',sColor);
%         hold on;
%     end
% end

end

function caPlotTrajectory(traForward)
 for cnt = 1:1:length(traForward)
        
         colorID = mod(cnt, 6);

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

        plot(round(traForward{cnt}(:,1)),round(traForward{cnt}(:,2)),'.-','LineWidth',4,'Color',sColor);
        hold on;
 end
end


function [orderNum] = caIndexOrderNum(intersectPtsIndex)
clear sortV;
clear sortI;
dataCol = intersectPtsIndex;
[m,n] = size(dataCol);
intersectPtsIndex2 = zeros(m,n);
for i = 1:1:m
    [sortV, sortI] =  sort(dataCol(i,:),'descend');
    numNonZero = length(find(sortV>0));
    sortOrder = zeros(1,numNonZero);
    
    for j = numNonZero:-1:2
        if j == numNonZero
            sortOrder(1,j) = 1;
        end
        if sortV(1,j-1) == sortV(1,j)%same index share the same order number
            sortOrder(1,j-1) = sortOrder(1,j);
        else
            sortOrder(1,j-1) = sortOrder(1,j) + 1;
        end   
    end
    sortV(1,1:numNonZero) = sortOrder;
    intersectPtsIndex2(i,sortI) = sortV;
end
orderNum = intersectPtsIndex2;

end



function [sortColNonZero] = caSortInterPts(intersectPtsTanIndex2)
dataCol = intersectPtsTanIndex2;

[iRow,iCol] = size(dataCol);
sortCol = zeros(5,2*iCol); %suppose there are at most 5 normal trajectories in one column.

indexVisited = zeros(1,iCol);
indexVisited = (sum(dataCol) == 0);

%%first round
cnt = 1;
while 1

    indexUnvisited = (indexVisited == 0);
    unVisitedMask = repmat(indexUnvisited,iRow,1);
    unVistedData = dataCol.*unVisitedMask;
    [subSort{cnt},indexVisited] = caFunSortPts(unVistedData);
    
    clear unVisitedCol;
     [~,unVisitedCol] = find(indexVisited == 0);
    if isempty(unVisitedCol)
        break;
    end
    
    cnt = cnt + 1;
    
end


%%%%% put the subsort together %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sortCol(1,(iCol - length(subSort{1}) + 1): iCol) = subSort{1};
if length(subSort) == 1
   sortColNonZero{1} = subSort{1};
end
if length(subSort) > 1
   %%%%%
   lengthNum = zeros(1,length(subSort));
   for i = 1:length(subSort)
       lengthNum(i) = length(subSort{i});
   end
   [~,iLengthDes] = sort(lengthNum(2:end),'descend');
   
   subSortDes{1} = subSort{1};
   for i = 1:length(iLengthDes)
       subSortDes{i+1} = subSort{iLengthDes(i)+1};
   end
   subSort = subSortDes;
   %%%%%

   indexVisited = zeros(1,length(subSort));
   indexVisited(1) = 1;
  
crossI = 0;
newNumUnvisited = 0;
oldNumUnvisited = 0;
while 1
   
    clear unVisitedCol;
    [~,unVisitedCol] = find(indexVisited == 0);
    if isempty(unVisitedCol)
       nonZeroCol = find(sum(sortCol)>0);
       nonZeroRow = find(sum(sortCol.')>0);
       sortNonZero = sortCol(nonZeroRow,nonZeroCol);
       sortColNonZero{1} = sortNonZero;
       break;
    end
    
    newNumUnvisited = length(unVisitedCol);
    if oldNumUnvisited == newNumUnvisited
       crossI = crossI + 1;
    else
       crossI = 0;
    end
    if crossI > 3
       nonZeroCol = find(sum(sortCol)>0);
       nonZeroRow = find(sum(sortCol.')>0);
       sortNonZero = sortCol(nonZeroRow,nonZeroCol);
       sortColNonZero{1} = sortNonZero;
        if ~isempty(unVisitedCol)
            for cnt = 1:1:length(unVisitedCol)
                sortColNonZero{cnt+1} = subSort{unVisitedCol(cnt)};
            end
        end
        break;
    end
    oldNumUnvisited = newNumUnvisited;
   
   for cnt = 1:length(unVisitedCol)
       
       i = unVisitedCol(cnt);
       if crossI < length(subSort{i})
           subSortIndex = subSort{i};
           dataColR = dataCol(:,subSortIndex(1+crossI));
           dataColL = dataCol(:,subSortIndex(end-crossI));
           
           iRNearestIndex = 0;
           iLNearestIndex = 0;
           
           validNumL = 0;
           validNumR = 0;
           
           nonZeroCol = find(sum(sortCol)>0);
           nonZeroRow = find(sum(sortCol.')>0);
           sortNonZero = sortCol(nonZeroRow,nonZeroCol);
           [rowN,colN] = size(sortNonZero);
           
           flagBreak1 = false;
           
           for cntCol = colN:-1:1
               if flagBreak1
                   flagBreak1 = false;
                   break;
               end
               
               for cntRow = rowN:-1:1
                   if sortNonZero(cntRow,cntCol) > 0
                       dataColStatic = dataCol(:,sortNonZero(cntRow,cntCol));

                       maskR = (dataColR > 0).*(dataColStatic>0);
                       diffR = maskR.*(dataColR - dataColStatic);

                       validNumR = length(find(diffR == 1)) + length(find(diffR + maskR == 1));

                       if (validNumR > 0.7*sum(maskR))
                          iRNearestIndex = sortNonZero(cntRow,cntCol);
                          flagBreak1 = true;
                          break;
                       end  
                   end
               end
           end
    
           for cntCol = 1:1:colN
               if flagBreak1
                   flagBreak1 = false;
                   break;
               end
               for cntRow = rowN:-1:1
                   if sortNonZero(cntRow,cntCol) > 0
                       dataColStatic = dataCol(:,sortNonZero(cntRow,cntCol));

                       maskL = (dataColL > 0).*(dataColStatic>0);
                       diffL = maskL.*( dataColStatic - dataColL);

                       validNumL = length(find(diffL == 1)) + length(find(diffL + maskL == 1));

                       if (validNumL > 0.7*sum(maskL))
                          iLNearestIndex = sortNonZero(cntRow,cntCol);
                          flagBreak1 = true;
                          break;
                       end  
                   end
               end
           end

           
           if (validNumR>= validNumL) && (iRNearestIndex > 0)    
               indexVisited(i) =1;
               subSortCut = subSort{i}(1,1+crossI:end);
               
              [indexR, indexCol] = find(sortCol==iRNearestIndex); 
               colBegin = indexCol + 1;
               colEnd = colBegin + length(subSortCut) - 1;
               
            if ~isempty(sortCol(sortCol(:,colBegin)>0, colBegin))
                 colIndex = sortCol(sortCol(:,colBegin)>0, colBegin);
                 dataColMask2 = dataCol(:,subSortCut(1))>0;
                 interFlag = zeros(1,length(colIndex));
                 for iCol = 1:1:length(colIndex)
                     dataColMask1 = dataCol(:,colIndex(iCol))>0;
                     interFlag(iCol) = sum(dataColMask1.*dataColMask2);
                 end
                 if sum(interFlag) == 0
                     if sum (sortCol(2,colBegin:colEnd)) == 0
                        sortCol(2,colBegin:colEnd) = subSortCut;
                     elseif sum(sortCol(3,colBegin:colEnd)) == 0
                        sortCol(3,colBegin:colEnd) = subSortCut;
                     elseif sum(sortCol(4,colBegin:colEnd)) == 0
                        sortCol(4,colBegin:colEnd) = subSortCut;
                     elseif sum(dataColMask1.*dataColMask2) == 0
                        sortCol(5,colBegin:colEnd) = subSortCut;  
                     end
                 end
            else
                 if sum (sortCol(2,colBegin:colEnd)) == 0
                    sortCol(2,colBegin:colEnd) = subSortCut;
                 elseif sum(sortCol(3,colBegin:colEnd)) == 0
                    sortCol(3,colBegin:colEnd) = subSortCut;
                 elseif sum(sortCol(4,colBegin:colEnd)) == 0
                    sortCol(4,colBegin:colEnd) = subSortCut;
                 elseif sum(dataColMask1.*dataColMask2) == 0
                    sortCol(5,colBegin:colEnd) = subSortCut;  
                 end    
            end
               
%              if sum (sortCol(2,colBegin:colEnd)) == 0
%                  if sortCol(1,colBegin)>0
%                      dataColMask1 = dataCol(:,sortCol(1,colBegin))>0;
%                      dataColMask2 = dataCol(:,subSortCut(1))>0;
%                      if sum(dataColMask1.*dataColMask2) == 0
%                         sortCol(2,colBegin:colEnd) = subSortCut;
%                      end 
%                  else
%                     sortCol(2,colBegin:colEnd) = subSortCut;
%                  end
%              elseif sum(sortCol(3,colBegin:colEnd)) == 0
%                 if sortCol(2,colBegin)>0
%                  dataColMask1 = dataCol(:,sortCol(2,colBegin))>0;
%                  dataColMask2 = dataCol(:,subSortCut(1))>0;
%                  if sum(dataColMask1.*dataColMask2) == 0
%                     sortCol(3,colBegin:colEnd) = subSortCut;
%                  end
%                 else
%                     sortCol(3,colBegin:colEnd) = subSortCut;
%                 end
%                 
%              elseif sum(sortCol(4,colBegin:colEnd)) == 0
%                if sortCol(3,colBegin)>0
%                  dataColMask1 = dataCol(:,sortCol(3,colBegin))>0;
%                  dataColMask2 = dataCol(:,subSortCut(1))>0;
%                  if sum(dataColMask1.*dataColMask2) == 0
%                     sortCol(4,colBegin:colEnd) = subSortCut;
%                  end
%                else
%                    sortCol(4,colBegin:colEnd) = subSortCut;
%                end
%              elseif sum(sortCol(5,colBegin:colEnd)) == 0
%                 if sortCol(4,colBegin)>0
%                  dataColMask1 = dataCol(:,sortCol(4,colBegin))>0;
%                  dataColMask2 = dataCol(:,subSortCut(1))>0;
%                  if sum(dataColMask1.*dataColMask2) == 0
%                     sortCol(5,colBegin:colEnd) = subSortCut;
%                  end
%                 else
%                    sortCol(5,colBegin:colEnd) = subSortCut; 
%                 end
%              end   
             
           end
           
           if (validNumR < validNumL) && (iLNearestIndex > 0)
               indexVisited(i) = 1;
               subSortCut = subSort{i}(1,1:end - crossI);
              [indexR, indexCol] = find(sortCol==iLNearestIndex); 
               colEnd = indexCol - 1;
               colBegin = colEnd - length(subSortCut)+1;
             if ~isempty(sortCol(sortCol(:,colEnd)>0, colEnd))
                 colIndex = sortCol(sortCol(:,colEnd)>0, colEnd);
                 dataColMask2 = dataCol(:,subSortCut(end))>0;
                 interFlag = zeros(1,length(colIndex));
                 for iCol = 1:1:length(colIndex)
                     dataColMask1 = dataCol(:,colIndex(iCol))>0;
                     interFlag(iCol) = sum(dataColMask1.*dataColMask2);
                 end
                 if sum(interFlag) == 0
                     if sum (sortCol(2,colBegin:colEnd)) == 0
                         sortCol(2,colBegin:colEnd) = subSortCut;
                     elseif sum(sortCol(3,colBegin:colEnd)) == 0
                         sortCol(3,colBegin:colEnd) = subSortCut;  
                     elseif sum(sortCol(4,colBegin:colEnd)) == 0
                         sortCol(4,colBegin:colEnd) = subSortCut; 
                     elseif sum(sortCol(5,colBegin:colEnd)) == 0
                         sortCol(5,colBegin:colEnd) = subSortCut;
                     end
                 end
            else
                 if sum (sortCol(2,colBegin:colEnd)) == 0
                     sortCol(2,colBegin:colEnd) = subSortCut;
                 elseif sum(sortCol(3,colBegin:colEnd)) == 0
                     sortCol(3,colBegin:colEnd) = subSortCut;  
                 elseif sum(sortCol(4,colBegin:colEnd)) == 0
                     sortCol(4,colBegin:colEnd) = subSortCut; 
                 elseif sum(sortCol(5,colBegin:colEnd)) == 0
                     sortCol(5,colBegin:colEnd) = subSortCut;
                 end
             end
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%              if sum (sortCol(2,colBegin:colEnd)) == 0
%                   if sortCol(1,colEnd)>0
%                      dataColMask1 = dataCol(:,sortCol(1,colEnd))>0;
%                      dataColMask2 = dataCol(:,subSortCut(end))>0;
%                      if sum(dataColMask1.*dataColMask2) == 0
%                         sortCol(2,colBegin:colEnd) = subSortCut;
%                      end 
%                  else
%                     sortCol(2,colBegin:colEnd) = subSortCut;
%                  end         
%              elseif sum(sortCol(3,colBegin:colEnd)) == 0
%                 if sortCol(2,colEnd)>0
%                  dataColMask1 = dataCol(:,sortCol(2,colEnd))>0;
%                  dataColMask2 = dataCol(:,subSortCut(end))>0;
%                  if sum(dataColMask1.*dataColMask2) == 0
%                     sortCol(3,colBegin:colEnd) = subSortCut;
%                  end
%                 else
%                     sortCol(3,colBegin:colEnd) = subSortCut;
%                 end
%                 
%              elseif sum(sortCol(4,colBegin:colEnd)) == 0
%                 if sortCol(3,colEnd)>0
%                  dataColMask1 = dataCol(:,sortCol(3,colEnd))>0;
%                  dataColMask2 = dataCol(:,subSortCut(end))>0;
%                  if sum(dataColMask1.*dataColMask2) == 0
%                     sortCol(4,colBegin:colEnd) = subSortCut;
%                  end
%                 else
%                     sortCol(4,colBegin:colEnd) = subSortCut;
%                 end
%                 
%              elseif sum(sortCol(5,colBegin:colEnd)) == 0
%                if sortCol(4,colEnd)>0
%                  dataColMask1 = dataCol(:,sortCol(4,colEnd))>0;
%                  dataColMask2 = dataCol(:,subSortCut(end))>0;
%                  if sum(dataColMask1.*dataColMask2) == 0
%                     sortCol(5,colBegin:colEnd) = subSortCut;
%                  end
%                else
%                     sortCol(5,colBegin:colEnd) = subSortCut;
%                end
%                 
%              end   
             
           end
                   
      end
   end
   
%     %%Debug Begin
%     nonZeroCol = find(sum(sortCol)>0);
%     nonZeroRow = find(sum(sortCol.')>0);
%     sortNonZero = sortCol(nonZeroRow,nonZeroCol);
%     
%     [rowNum,colNum] = size(sortNonZero);
%     [dataRowNum,~] = size(dataCol);
%     sortIndexTemp = zeros(dataRowNum,colNum);
%     
%     for iCol = 1:colNum
%         for iRow = 1:rowNum
%             if sortNonZero(iRow,iCol)>0
%                 
%             sortIndexTemp(:,iCol) = sortIndexTemp(:,iCol) + dataCol(:,sortNonZero(iRow,iCol));
%             end
%         end   
%     end
%     figure(10);
%     imagesc(sortIndexTemp);
%     %%Debu End
    
end
end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% sortColNonZero{1} = sortNonZero;

end

function [sortColNonZero, indexVisited] = caFunSortPts(intersectPtsTanIndex2)

dataCol = intersectPtsTanIndex2;

[~,iCol] = size(dataCol);
sortCol = zeros(1,2*iCol); %suppose there are at most 5 normal trajectories in one column.

%find the longest col as the inital col
[~,iLongestCol]  = max(sum(dataCol > 0));

indexVisited = zeros(1,iCol);
indexVisited = (sum(dataCol) == 0);

iL = iCol + 1;
iR = iCol + 1;
%%first round
while 1
    clear unVisitedCol;
    [~,unVisitedCol] = find(indexVisited == 0);
    if isempty(unVisitedCol)
        break;
    end
    
    if iL == iR
       sortCol(1,iR) = iLongestCol;
       indexVisited(sortCol(1,iR)) = 1;
       clear unVisitedCol;
       [~,unVisitedCol] = find(indexVisited == 0);
    end
    
    validNumRMax = 0;
    validNumLMax = 0;
    
    iRNearestCol = 0;
    iLNearestCol = 0;
   
    for i = 1:1:length(unVisitedCol)
        
        maskR = (dataCol(:,unVisitedCol(i)) > 0).*(dataCol(:,sortCol(1,iR))>0);
        diffR = maskR.*(dataCol(:,unVisitedCol(i)) - dataCol(:,sortCol(1,iR)));
        
        maskL = (dataCol(:,unVisitedCol(i)) > 0).*(dataCol(:,sortCol(1,iL))>0);
        diffL = maskL.*(dataCol(:,sortCol(1,iL)) - dataCol(:,unVisitedCol(i)));
        
        validNumR = length(find(diffR == 1)) + length(find(diffR + maskR == 1));
        validNumL = length(find(diffL == 1)) + length(find(diffL + maskL == 1));
        
        if (validNumR > 0.7*sum(maskR)) && (validNumR > validNumRMax) && (validNumR >= validNumL)
            validNumRMax = validNumR;
            iRNearestCol = unVisitedCol(i);
            %iRNearestCol{length(iRNearestCol)+1} = unVisitedCol(i);

        end
        
        if (validNumL > 0.7*sum(maskL)) && (validNumL > validNumLMax) && (validNumL > validNumR)
            validNumLMax = validNumL; 
            iLNearestCol = unVisitedCol(i);
%             iLNearestCol{length(iLNearestCol)+1} = unVisitedCol(i);
        end
    end
    
    if validNumRMax > 0
        iR = iR + 1;
        sortCol(1,iR) = iRNearestCol;
         indexVisited(sortCol(1,iR)) = 1;
    end
    
    if validNumLMax > 0
       iL = iL - 1;
       sortCol(1,iL) = iLNearestCol;
        indexVisited(sortCol(1,iL)) = 1;
    end
 
    
%     sortColNonZero = sortCol(find(sortCol>0)); 
%     %%Debug Begin
%     sortColNonZeroTemp = sortCol(find(sortCol>0));    
%     sortIndexTemp = dataCol(:,sortColNonZeroTemp);
%     figure(10);
%     imagesc(sortIndexTemp);
%     %%Debu End
    
    if validNumRMax == 0 && validNumLMax == 0
        break;
    end
    
end

 if iR - (iCol+1) > (iCol+1) - iL
     halfSortCol = [sortCol(1,1:iCol),zeros(1,iCol)];
    if length(find(sortCol(halfSortCol>0)))>0
       indexVisited(sortCol(halfSortCol>0)) = 0;
    end
    sortCol(1,1:iCol) = zeros(1,iCol);
 elseif iR - (iCol+1) < (iCol+1) - iL
     halfSortCol = [zeros(1,iCol+1),sortCol(1,iCol+2:end)];
    if length(find(sortCol(halfSortCol>0)))>0
       indexVisited(sortCol(halfSortCol>0)) = 0;
    end
    sortCol(1,iCol+2:end) = zeros(1,iCol-1);
 end
 
   sortColNonZero = sortCol(find(sortCol>0)); 
    
%     %%Debug Begin
%     sortColNonZeroTemp = sortCol(find(sortCol>0));    
%     sortIndexTemp = dataCol(:,sortColNonZeroTemp);
%     figure(10);
%     imagesc(sortIndexTemp);
%     %%Debu End
    
end


%%%%%%%%%%%%%%% - old version of caSortInterPts - %%%%%%%%%%%%%%%%%%%%%%
% function [sortColNonZero] = caSortInterPts(intersectPtsTanIndex2)
% dataCol = intersectPtsTanIndex2;
% 
% [~,iCol] = size(intersectPtsTanIndex2);
% sortCol = zeros(5,2*iCol); %suppose there are at most 5 normal trajectories in one column.
% 
% indexVisited = zeros(1,iCol);
% indexVisited = (sum(intersectPtsTanIndex2) == 0);
% 
% iL = iCol + 1;
% iR = iCol + 1;
% %%first round
% while 1
%     clear unVisitedCol;
%     [~,unVisitedCol] = find(indexVisited == 0);
%     if isempty(unVisitedCol)
%         break;
%     end
%     
%     if iL == iR
%        sortCol(1,iR) = unVisitedCol(1);
%        indexVisited(sortCol(1,iR)) = 1;
%        clear unVisitedCol;
%        [~,unVisitedCol] = find(indexVisited == 0);
%     end
%     
%     validNumRMax = 0;
%     validNumLMax = 0;
%     iRNearestCol = 0;
%     iLNearestCol = 0;
%     for i = 1:1:length(unVisitedCol)
%         
%         maskR = (dataCol(:,unVisitedCol(i)) > 0).*(dataCol(:,sortCol(1,iR))>0);
%         diffR = maskR.*(dataCol(:,unVisitedCol(i)) - dataCol(:,sortCol(1,iR)));
%         
%         maskL = (dataCol(:,unVisitedCol(i)) > 0).*(dataCol(:,sortCol(1,iL))>0);
%         diffL = maskL.*(dataCol(:,sortCol(1,iL)) - dataCol(:,unVisitedCol(i)));
%         
%         validNumR = length(find(diffR == 1)) + length(find(diffR + maskR == 1));
%         validNumL = length(find(diffL == 1)) + length(find(diffL + maskL == 1));
%         
%         if (validNumR > 0.7*sum(maskR)) && (validNumR > validNumRMax) && (validNumR >= validNumL)
%             validNumRMax = validNumR;
%             iRNearestCol = unVisitedCol(i);
%         end
%         
%         if (validNumL > 0.7*sum(maskL)) && (validNumL > validNumLMax) && (validNumL > validNumR)
%             validNumLMax = validNumL; 
%             iLNearestCol = unVisitedCol(i);
%         end
%     end
%     
%     if validNumRMax > 0
%         iR = iR + 1;
%         sortCol(1,iR) = iRNearestCol;
%         indexVisited(sortCol(1,iR)) = 1;
%     end
%     
%     if validNumLMax > 0
%        iL = iL - 1;
%        sortCol(1,iL) = iLNearestCol;
%        indexVisited(sortCol(1,iL)) = 1;
%     end
%  
%     if validNumRMax == 0 && validNumLMax == 0
%         break;
%     end
%     
%     %%Debug Begin
%     sortColNonZeroTemp = sortCol(find(sortCol>0));    
%     sortIndexTemp = dataCol(:,sortColNonZeroTemp);
%     figure(10);
%     imagesc(sortIndexTemp);
%     %%Debu End
%     
% end
% 
% %second round, insert the unvisited col to the sorted col.
% flagStaticCol = 0;
% while 1
%     clear unVisitedCol;
%     clear visitedCol;
%     [~,unVisitedCol] = find(indexVisited == 0);
%     if isempty(unVisitedCol)
%         break;
%     end
% %if ~isempty(unVisitedCol)
%     for i = 1:1:length(unVisitedCol)
%         validNumRMax = 0;
%         validNumLMax = 0;
%         iRNearestIndex = 0;
%         iLNearestIndex = 0;
%         
%         visitedCol = sortCol(sortCol>0);
%         for j = 1:1:length(visitedCol)
%            maskR = (dataCol(:,unVisitedCol(i)) > 0).*(dataCol(:,visitedCol(j))>0);
%            diffR = maskR.*(dataCol(:,unVisitedCol(i)) - dataCol(:,visitedCol(j)));
%            
%            maskL = (dataCol(:,unVisitedCol(i)) > 0).*(dataCol(:,visitedCol(j))>0);
%            diffL = -diffR;
% 
%            validNumR = length(find(diffR == 1)) + length(find(diffR + maskR == 1));
%            validNumL = length(find(diffL == 1)) + length(find(diffL + maskL == 1));
%            
%            if (validNumR > 0.7*sum(maskR)) && (validNumR >= validNumRMax)
%               validNumRMax = validNumR;
%               iRNearestIndex = visitedCol(j);
%            end  
%            
%            if (validNumL > 0.7*sum(maskL)) && (validNumL > validNumLMax)
%               validNumLMax = validNumL;
%               iLNearestIndex = visitedCol(j);
%            end  
%            
%         end
% 
%         if (validNumRMax > 0) && (validNumRMax >= validNumLMax)
%            [iRow, iCol] = find(sortCol==iRNearestIndex); 
%            if iRow ==1
%                if sortCol(iRow,iCol+1)>0
%                maskR1 = dataCol(:,sortCol(iRow,iCol+1))>0;
%                maskR2 = dataCol(:,unVisitedCol(i))>0;
%                if ~isempty(find(maskR1.*maskR2 ==1))
%                   sortCol(:,iCol+2:end) = sortCol(:,iCol+1:end-1);
%                   sortCol(:,iCol+1) = zeros(5,1);
%                   sortCol(iRow,iCol+1) = unVisitedCol(i); 
%                else
%                   sortCol(iRow+1,iCol+1) = unVisitedCol(i);
%                   flagStaticCol = iCol;
%                end
%                else
%                   sortCol(:,iCol+2:end) = sortCol(:,iCol+1:end-1);
%                   sortCol(:,iCol+1) = zeros(5,1);
%                   sortCol(iRow,iCol+1) = unVisitedCol(i); 
%                end
%            else 
%                if iCol > flagStaticCol 
%                    if sortCol(iRow, iCol+1) == 0
%                       sortCol(iRow,iCol+1) = unVisitedCol(i);
%                    else
%                       sortCol(iRow,iCol+2:end) = sortCol(iRow,iCol+1:end-1);
%                       sortCol(iRow,iCol+1) = unVisitedCol(i);  
%                    end
%                else
%                    if sortCol(iRow,iCol-1)==0
%                       sortCol(iRow,iCol-1) = unVisitedCol(i);
%                    else
%                       sortCol(iRow,1:iCol-1) = sortCol(iRow,2:iCol);
%                       sortCol(iRow,iCol) = unVisitedCol(i); 
%                    end              
%                end
%            end
%                
%            indexVisited(unVisitedCol(i)) = 1;
%            
%         elseif (validNumLMax > 0) && (validNumLMax > validNumRMax)
%            [iRow, iCol] = find(sortCol==iLNearestIndex); 
%            if iRow ==1
%                if sortCol(iRow,iCol-1)>0
%                    maskL1 = dataCol(:,sortCol(iRow,iCol-1))>0;
%                    maskL2 = dataCol(:,unVisitedCol(i))>0;
%                    if ~isempty(find(maskL1.*maskL2 ==1))
%                       sortCol(:,iCol+1:end) = sortCol(:,iCol:end-1);
%                       sortCol(:,iCol) = zeros(5,1);
%                       sortCol(iRow,iCol) = unVisitedCol(i);
%                    else
%                       sortCol(iRow+1,iCol-1) = unVisitedCol(i);
%                        flagStaticCol = iCol;
%                    end
%                else
%                    sortCol(:,iCol+1:end) = sortCol(:,iCol:end-1);
%                    sortCol(:,iCol) = zeros(5,1);
%                    sortCol(iRow,iCol) = unVisitedCol(i);
%                end
%            else
%                if iCol > flagStaticCol
%                    if sortCol(iRow,iCol-1)==0
%                       sortCol(iRow,iCol-1) = unVisitedCol(i);
%                    else
%                        sortCol(iRow,iCol+1:end) = sortCol(iRow,iCol:end-1);
%                        sortCol(iRow,iCol) = unVisitedCol(i); 
%                    end       
%                else
%                    if sortCol(iRow,iCol-1)==0
%                       sortCol(iRow,iCol-1) = unVisitedCol(i);
%                    else
%                       sortCol(iRow,1:iCol-1) = sortCol(iRow,2:iCol);
%                       sortCol(iRow,iCol) = unVisitedCol(i); 
%                    end       
%                end
%            end
%            indexVisited(unVisitedCol(i)) = 1;
%         end
%         
%     %%Debug Begin
%     nonZeroCol = find(sum(sortCol)>0);
%     nonZeroRow = find(sum(sortCol.')>0);
%     sortNonZero = sortCol(nonZeroRow,nonZeroCol);
%     
%     [rowNum,colNum] = size(sortNonZero);
%     [dataRowNum,~] = size(dataCol);
%     sortIndexTemp = zeros(dataRowNum,colNum);
%     
%     for iCol = 1:colNum
%         for iRow = 1:rowNum
%             if sortNonZero(iRow,iCol)>0
%             sortIndexTemp(:,iCol) = sortIndexTemp(:,iCol) + dataCol(:,sortNonZero(iRow,iCol));
%             end
%         end   
%     end
% %     sortColNonZeroTemp = sortCol(find(sortCol>0));    
% %     sortIndexTemp = dataCol(:,sortColNonZeroTemp);
%     figure(10);
%     imagesc(sortIndexTemp);
%     %%Debu End
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  sortColNonZero = sortCol(find(sortCol>0));
%  
%  end

function [intersectPts,intersectPtsTanIndex,intersectPtsNorIndex] = caComputIntersetPts(trajectoryTan,trajectoryNor)

numTan = length(trajectoryTan);
numNor = length(trajectoryNor);

intersectPts = zeros(numTan, numNor, 2);
intersectPtsTanIndex = zeros(numTan, numNor);
intersectPtsNorIndex = zeros(numTan, numNor);

for iRow = 1:1:numTan 
for iCol = 1:1:numNor
    [PtsTemp,ia,ib] = intersect(trajectoryTan{iRow},trajectoryNor{iCol},'rows');
    if ~isempty(PtsTemp)
        if length(PtsTemp(:,1)) == 1
           intersectPts(iRow,iCol,1) = PtsTemp(1,1);
           intersectPts(iRow,iCol,2) = PtsTemp(1,2);
           intersectPtsTanIndex(iRow,iCol) = ia(1,1);
           intersectPtsNorIndex(iRow,iCol) = ib(1,1);
        elseif length(PtsTemp(:,1)) == 2
           intersectPts(iRow,iCol,1) = 0.5 * (PtsTemp(1,1) + PtsTemp(2,1));
           intersectPts(iRow,iCol,2) = 0.5 * (PtsTemp(1,2) + PtsTemp(2,2));
           intersectPtsTanIndex(iRow,iCol) = 0.5 * (ia(1,1) + ia(2,1));
           intersectPtsNorIndex(iRow,iCol) = 0.5 * (ib(1,1) + ib(2,1));
        end
    else %相交位置为0.5,0.5
       clear trajectoryTemp;
       trajectoryTemp(:,1) = trajectoryTan{iRow}(:,1);
       trajectoryTemp(:,1) = trajectoryTemp(:,1) + 1;
       trajectoryTemp(:,2) = trajectoryTan{iRow}(:,2);
       [PtsTemp1,ia1,ib1] = intersect(trajectoryTemp,trajectoryNor{iCol},'rows'); 
       clear trajectoryTemp;
       trajectoryTemp(:,1) = trajectoryTan{iRow}(:,1);
       trajectoryTemp(:,1) =  trajectoryTemp(:,1) - 1;
       trajectoryTemp(:,2) = trajectoryTan{iRow}(:,2);
       [PtsTemp2,ia2,ib2] = intersect(trajectoryTemp,trajectoryNor{iCol},'rows'); 
       if ~isempty(PtsTemp1) && ~isempty(PtsTemp2)
           if length(PtsTemp1(:,1)) == 2
               intersectPts(iRow,iCol,1) = 0.5*(PtsTemp1(1,1) + PtsTemp1(2,1));
               intersectPts(iRow,iCol,2) = 0.5*(PtsTemp1(1,2) + PtsTemp1(2,2));
               intersectPtsTanIndex(iRow,iCol) = 0.5*(ia1(1,1) + ia1(2,1));
               intersectPtsNorIndex(iRow,iCol) = 0.5*(ib1(1,1) + ib1(2,1));
           elseif length(PtsTemp2(:,2)) == 2
              intersectPts(iRow,iCol,1) = 0.5*(PtsTemp2(1,1) + PtsTemp2(2,1));
              intersectPts(iRow,iCol,2) = 0.5*(PtsTemp2(1,2) + PtsTemp2(2,2));
              intersectPtsTanIndex(iRow,iCol) = 0.5*(ia2(1,1) + ia2(2,1));
              intersectPtsNorIndex(iRow,iCol) = 0.5*(ib2(1,1) + ib2(2,1));
           elseif length(PtsTemp1(:,1))==1 && length(PtsTemp2(:,1))==1 
              intersectPts(iRow,iCol,1) = 0.5*(PtsTemp1(1,1) + PtsTemp2(1,1));
              intersectPts(iRow,iCol,2) = 0.5*(PtsTemp1(1,2) + PtsTemp2(1,2));
              intersectPtsTanIndex(iRow,iCol) = 0.5*(ia1 + ia2);
              intersectPtsNorIndex(iRow,iCol) = 0.5*(ib1 + ib2);
           end
       end
    end
    
end
end
end


function [curlMap,divMap,curlMat,divMat,curlFeatureVector,divFeatureVector] = caCDdescriptor(uv,sortIntersectPts)
clear magnitudeFlow;
magnitudeFlow(:,:,1) = abs(uv(:,:,1)) + abs(uv(:,:,2));
magnitudeFlow(:,:,2) = abs(uv(:,:,1)) + abs(uv(:,:,2));
IndexMask = (magnitudeFlow(:,:,1)~= 0);
se = strel('disk',1);    
erodedIndexMask = imerode(IndexMask,se);

uvNorm = caFlowNormlize(uv);
u = uvNorm(:,:,1);
v = uvNorm(:,:,2);

[~,caz] = curl (u,v);
cazROI = erodedIndexMask.*caz;
%  cazYFlip = flipud(cazROI);
curlMap = cazROI;

div = divergence(u,v);
divROI = erodedIndexMask.*div;
divMap = divROI;

%curlMatrix & divMatrix
[matR,matC,~] = size(sortIntersectPts);
curlMat = zeros(matR,matC);
divMat = zeros(matR,matC);

[m,n,z] = size(uv);
[Xmesh Ymesh] = meshgrid(1:n, 1:m);

curlMat =  interp2(Xmesh,Ymesh,cazROI,sortIntersectPts(:,:,1), sortIntersectPts(:,:,2),'nearest',0);
divMat =  interp2(Xmesh,Ymesh,divROI,sortIntersectPts(:,:,1), sortIntersectPts(:,:,2),'nearest',0);

curlFeatureVector = sum(curlMat.');
divFeatureVector = sum(divMat);
end


function [curlMap,divMap,curlInteDescriptor, divInteDescriptorAll,divInteDescriptorShort,curlInteMap,divInteMap, maskCurlMap,maskDivMap] = caCDIntegrationDescriptor(uv,sortRowNonZero,sortColNonZero,sortColNonZeroCut,trajectoryTan,trajectoryNor,intersectPts)
clear magnitudeFlow;
magnitudeFlow(:,:,1) = abs(uv(:,:,1)) + abs(uv(:,:,2));
magnitudeFlow(:,:,2) = abs(uv(:,:,1)) + abs(uv(:,:,2));
IndexMask = (magnitudeFlow(:,:,1)~= 0);
se = strel('disk',1);    
erodedIndexMask = imerode(IndexMask,se);

uvNorm = caFlowNormlize(uv);
u = uvNorm(:,:,1);
v = uvNorm(:,:,2);

[caz,~] = curl (u,v);
cazROI = erodedIndexMask.*caz;
%  cazYFlip = flipud(cazROI);
curlMap = cazROI;

div = divergence(u,v);
divROI = erodedIndexMask.*div;
divMap = divROI;

%CD integration descriptor
[m,n,z] = size(uv);
[Xmesh Ymesh] = meshgrid(1:n, 1:m);

[numTralets,lengthC] = size(sortRowNonZero);
if lengthC > 0
curlInteDescriptor = zeros(1,lengthC);
for i = 1:1:lengthC
    for j = 1:1:numTralets 
        traNum = sortRowNonZero(j,i);
        if traNum > 0
            curlOnTra = interp2(Xmesh,Ymesh,cazROI,trajectoryTan{traNum}(:,1), trajectoryTan{traNum}(:,2),'nearest',0); 
            curlInteDescriptor(1,i) = curlInteDescriptor(1,i) + sum(curlOnTra);
        end
%         curlOnTra = interp2(Xmesh,Ymesh,cazROI,trajectoryTan{traNum}(:,1), trajectoryTan{traNum}(:,2),'nearest',0);
%         curlInteDescriptor(1,i) = sum(curlOnTra);
    end
end
end



% [numTralets,lengthD] = size(sortColNonZero);
[numTralets,lengthD] = size(sortColNonZeroCut);
if lengthD > 0
divInteDescriptorShort = zeros(1,lengthD);
for i = 1:1:lengthD
    for j = 1:1:numTralets
        traNum = sortColNonZeroCut(j,i);
        if traNum>0
            divOnTra = interp2(Xmesh,Ymesh,divROI,trajectoryNor{traNum}(:,1), trajectoryNor{traNum}(:,2),'nearest',0);
            divInteDescriptorShort(1,i) = divInteDescriptorShort(1,i) + sum(divOnTra);
        end
    end
end
end


[numTralets,lengthD] = size(sortColNonZero);
%[numTralets,lengthD] = size(sortColNonZeroCut);
if lengthD > 0
divInteDescriptorAll = zeros(1,lengthD);
for i = 1:1:lengthD
    for j = 1:1:numTralets
        traNum = sortColNonZero(j,i);
        if traNum>0
            divOnTra = interp2(Xmesh,Ymesh,divROI,trajectoryNor{traNum}(:,1), trajectoryNor{traNum}(:,2),'nearest',0);
            divInteDescriptorAll(1,i) = divInteDescriptorAll(1,i) + sum(divOnTra);
        end
    end
end
end

%CD integration Map

[numTraTanlets,iRow] = size(sortRowNonZero);
[numTraNorlets,iCol] = size(sortColNonZero);

if iRow>0 && iCol>0
    
[~,intPtsN,~] = size(intersectPts);

sortIntersectPts = zeros(iRow,iCol,2);
sortIntersectPtsRow = zeros(iRow,intPtsN,2);

for i = 1:1:iRow
    for j=1:1:numTraTanlets
        if sortRowNonZero(j,i)>0
           sortIntersectPtsRow(i,:,:) = sortIntersectPtsRow(i,:,:) + intersectPts(sortRowNonZero(j,i),:,:);
        end
    end
end

for i = 1:1:iCol
    for j = 1:1:numTraNorlets
         if sortColNonZero(j,i)>0
           sortIntersectPts(:,i,:) = sortIntersectPts(:,i,:) + sortIntersectPtsRow(:,sortColNonZero(j,i),:);
        end
    end
end

% figure;
% imagesc(sortIntersectPts(:,:,1));
% figure;
% imagesc(sortIntersectPts(:,:,2));
% 
% 
%  sortIntersectPts = intersectPts(sortRowNonZero,sortColNonZero,:);
    
[m,n,z] = size(uv);
[Xmesh Ymesh] = meshgrid(1:n, 1:m);

curlInteMap = zeros(iRow,iCol);
divInteMap = zeros(iRow,iCol);

maskCurlMap = -1 * ones(iRow,iCol);
maskDivMap = -1 * ones(iRow,iCol);

for i = 1:1:iRow
    for j = 1:1:iCol
      
        Pts(1,:) = sortIntersectPts(i,j,:);
        if Pts(1,1)>0 || Pts(1,2)>0
           Pts(2,:) = ceil(Pts(1,:));
           Pts(3,:) = floor(Pts(1,:));
           Pts(4,1) = ceil(Pts(1,1)); Pts(4,2) = floor(Pts(1,2));
           Pts(5,1) = floor(Pts(1,1)); Pts(5,2) = ceil(Pts(1,2));
           Pts(6,1) = Pts(2,1) - 1;Pts(6,2) = Pts(2,2);
           Pts(7,1) = Pts(2,1) + 1;Pts(7,2) = Pts(2,2);
           Pts(8,1) = Pts(3,1) - 1;Pts(8,2) = Pts(3,2);
           Pts(9,1) = Pts(3,1) + 1;Pts(9,2) = Pts(3,2);
           
           for iTraTanlets = 1:1:numTraTanlets
               if sortRowNonZero(iTraTanlets,i)>0
                  traTanNum = sortRowNonZero(iTraTanlets,i);
                  [~,iTan,~] = intersect(trajectoryTan{traTanNum},Pts,'rows'); 
                  if ~isempty(iTan)
                      curlOnTra = interp2(Xmesh,Ymesh,cazROI,trajectoryTan{traTanNum}(:,1), trajectoryTan{traTanNum}(:,2),'nearest',0);
                      
                      maskCurl = zeros(1,length(trajectoryTan{traTanNum}));
                      maskCurl(1,1:max(iTan)) = 1;
                      curlInteMap(i,j) = maskCurl * curlOnTra;
                      maskCurlMap(i,j) = traTanNum;
                      break;
                  end
               end
           end
           
           for iTraNorlets = 1:1:numTraNorlets
               if sortColNonZero(iTraNorlets,j)>0;
                 traNorNum = sortColNonZero(iTraNorlets,j);
                 [~,iNor,~] = intersect(trajectoryNor{traNorNum},Pts,'rows');   
                 if ~isempty(iNor)
                      divOnTra = interp2(Xmesh,Ymesh,divROI,trajectoryNor{traNorNum}(:,1), trajectoryNor{traNorNum}(:,2),'nearest',0);
       
                     maskDiv = zeros(1,length(trajectoryNor{traNorNum}));
                     maskDiv(1,1:max(iNor)) = 1;
                     divInteMap(i,j) = maskDiv * divOnTra;
                     maskDivMap(i,j) = traNorNum;
                      break;
                  end
               end 
           end
        end
        
    end

end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sortColNonZero,sortRowNonZero] = caSortRowColInterPts(intersectPtsTanIndex2,intersectPtsNorIndex2)
[sortColTemp] = caSortInterPts(intersectPtsTanIndex2);
sortRow = {};
sortCol = {};
for cnt = 1:1:length(sortColTemp)
sortNonZero = sortColTemp{cnt};
[~,sizeN] = size(sortNonZero);
    if sizeN > 10  %ignore small part
        indexNonZero = sortNonZero(sortNonZero>0);
        dataCol = intersectPtsNorIndex2;
        [~,colNum] = size(dataCol);
        indexZero = zeros(1,colNum);
        indexZero(1,indexNonZero) = 1;
        indexColZero = find(indexZero==0);
        dataCol(:,indexColZero) = 0;
        [sortRowTemp] = caSortInterPts(dataCol.');
        [~,sizeN] = size(sortRowTemp{1});
        if sizeN > 2
            sortRow{length(sortRow)+1} = sortRowTemp{1};
            sortCol{length(sortCol)+1} = sortColTemp{cnt};
        end
    end
end
sortColNonZero = sortCol;
sortRowNonZero = sortRow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function caShowSortIndex(dataCol,sortColNonZero)
   for cnt = 1:1:length(sortColNonZero)
     sortNonZero = sortColNonZero{cnt};
    [rowNum,colNum] = size(sortNonZero);
    [dataRowNum,~] = size(dataCol);
    sortIndexTemp = zeros(dataRowNum,colNum);
    
    for iCol = 1:colNum
        for iRow = 1:rowNum
            if sortNonZero(iRow,iCol)>0
            sortIndexTemp(:,iCol) = sortIndexTemp(:,iCol) + dataCol(:,sortNonZero(iRow,iCol));
            end
        end   
    end
    figure();
    imagesc(sortIndexTemp);
    title('sort index results');
   end
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function  [sortColNonZeroCut] = caSortColIndexCut(dataCol,sortColNonZero)
   sortColNonZeroCut = {};
   for cnt = 1:1:length(sortColNonZero)
     sortNonZero = sortColNonZero{cnt};
     MaskNonZero = sortNonZero>0;
     sumMask = sum(MaskNonZero);
     MaskNonOverlap = (sumMask ==1);
     
     lengthMat = zeros(1,length(MaskNonOverlap));
     for i = 1:1:length(MaskNonOverlap)
         if MaskNonOverlap(i)>0
            colIndex = sum(sortNonZero(:,i));
            lengthMat(i) = sum(dataCol(:,colIndex)>0);
         end
     end
     
     [lenMax,Imax] = max(lengthMat);
     
     IL = Imax;
     IR = Imax;
     
     while IL > 1
         if MaskNonOverlap(IL-1) == 1 && lengthMat(IL-1) > 0.9*lenMax
             IL = IL - 1;
         else
            break;
         end
     end
     
     while IR < length(MaskNonOverlap) && lengthMat(IR+1) > 0.9*lenMax
          if MaskNonOverlap(IR+1) == 1
             IR = IR + 1;
         else
            break;
         end
     end
     
    sortColNonZeroCut{cnt} =  sortNonZero(:,IL:IR);
   end
 
 end

