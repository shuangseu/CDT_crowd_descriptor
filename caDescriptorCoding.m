function [motionCode] = caDescriptorCoding(curlDescriptor,divDescriptor,codeLength,type)

disp('Motion coding...');

%codeLength = 6;
%scaleLength = 140;
scaleLength = (codeLength+1)*20;

curlDesScale = imresize(curlDescriptor,[1,scaleLength],'nearest');
divDesScale = imresize(divDescriptor,[1,scaleLength],'nearest');

curlDesScale = medfilt1(curlDesScale,5);
divDesScale = medfilt1(divDesScale,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure;
%    bar(curlDesScale);
%    title('curl descriptor scale');
%    
%    figure;
%    bar(divDesScale);
%    title('divergence descriptor scale');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pooling, sliding window with 50% overlap
if type == 1
    curlCode = zeros(1,codeLength);
    divCode = zeros(1,codeLength);
    for cnt = 1:1:codeLength
        iBegin = (cnt - 1)*20 + 1;
        iEnd = iBegin + 39;

        dataCurl = curlDesScale(1,iBegin:iEnd);
        maskPositive = dataCurl >= 0;
        maskNegtive = dataCurl <= 0;
        if sum(maskPositive) >= sum(maskNegtive)
           curlCode(1,cnt) = max(dataCurl);
        else
           curlCode(1,cnt) = min(dataCurl); 
        end

        dataDiv = divDesScale(1,iBegin:iEnd);
        maskPositive = dataDiv >= 0;
        maskNegtive = dataDiv <= 0;
        if sum(maskPositive) >= sum(maskNegtive)
           divCode(1,cnt) = max(dataDiv);
        else
           divCode(1,cnt) = min(dataDiv); 
        end
    end

    motionCode = [curlCode,divCode];
end

%Pooling, sliding window with no overlap
if type == 2
 
    scaleLength = (codeLength)*20;

    curlDesScale = imresize(curlDescriptor,[1,scaleLength],'nearest');
    divDesScale = imresize(divDescriptor,[1,scaleLength],'nearest');

    curlDesScale = medfilt1(curlDesScale,5);
    divDesScale = medfilt1(divDesScale,5);

    curlCode = zeros(1,codeLength);
    divCode = zeros(1,codeLength);
    for cnt = 1:1:codeLength
        iBegin = (cnt - 1)*20 + 1;
        iEnd = iBegin + 19;

        dataCurl = curlDesScale(1,iBegin:iEnd);
        maskPositive = dataCurl >= 0;
        maskNegtive = dataCurl <= 0;
        if sum(maskPositive) >= sum(maskNegtive)
           curlCode(1,cnt) = max(dataCurl);
        else
           curlCode(1,cnt) = min(dataCurl); 
        end

        dataDiv = divDesScale(1,iBegin:iEnd);
        maskPositive = dataDiv >= 0;
        maskNegtive = dataDiv <= 0;
        if sum(maskPositive) >= sum(maskNegtive)
           divCode(1,cnt) = max(dataDiv);
        else
           divCode(1,cnt) = min(dataDiv); 
        end
    end
    motionCode = [curlCode,divCode];
end

%just down-sampling
if type == 3
    curlCode = imresize(curlDescriptor,[1,codeLength],'nearest');
    divCode = imresize(divDescriptor,[1,codeLength],'nearest');
    motionCode = [curlCode,divCode];
end

%average pooling
if type == 4
    scaleLength = (codeLength)*20;

    curlDesScale = imresize(curlDescriptor,[1,scaleLength],'nearest');
    divDesScale = imresize(divDescriptor,[1,scaleLength],'nearest');

    curlDesScale = medfilt1(curlDesScale,5);
    divDesScale = medfilt1(divDesScale,5);

    curlCode = zeros(1,codeLength);
    divCode = zeros(1,codeLength);
    for cnt = 1:1:codeLength
        iBegin = (cnt - 1)*20 + 1;
        iEnd = iBegin + 19;

        dataCurl = curlDesScale(1,iBegin:iEnd);
        curlCode(1,cnt) = mean(dataCurl);

        dataDiv = divDesScale(1,iBegin:iEnd);
        divCode(1,cnt) = mean(dataDiv);
    end
    motionCode = [curlCode,divCode];
end


% for cnt = 1:1:length(motionCode)
%    figure;
%    bar(motionCode{cnt});
%    title('motion code');
% end

end