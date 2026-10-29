function [pointAdjX,labelNum] = getPointDesV2(wholeMatch,K,nameta,rotationAngSet,nameta1,theta)
                                              
%%
Knum = length(K);
num = size(wholeMatch,1);
X = wholeMatch(:,1:2);
Y = wholeMatch(:,4:5);
D1 = pdist2(X,X);
D2 = pdist2(Y,Y);
[D1s,I1] = sort(D1,2,'ascend');
radi = mean(D1s(:))*0.6;
D1Op=cell(Knum,1);I1op=cell(Knum,1);
for i=1:Knum
    D1Op{i} = D1s(:,2:K(i)+1);
    fitDist1=zeros(size(D1Op{i}));
    fitDist1(find(D1Op{i}<radi))=1;
    D1Op{i} = D1Op{i}.*fitDist1;
    I1op{i} = I1(:,2:K(i)+1);
end
%%
pointAdjX = cell(num,1);labelNum=1;
% initial
for i=1:num
    pointAdjX{i}.vecCluLabel=[];% cluster index number 
    pointAdjX{i}.adjPoi=[];% neighbor feature point index
    pointAdjX{i}.simiVec=[];% similar 相似特征匹配的点索引
end
for i=1:num
    for j=1:size(rotationAngSet,2)
        %tmpPoiW the neighborhood similarity for current point i at current angle
        %tmpPoiNum the multi-scale statistical count within the neighborhood of current point i at current angle j
        tmpPoiW = zeros(num,1);tmpPoiNum=zeros(num,1);
        for kn=1:Knum
            curAdjIndx = I1op{kn}(i,find(D1Op{kn}(i,:)~=0 & I1op{kn}(i,:)));
            %         curYAdjIndex = I2op(i,find(D2Op(i,:)~=0 & I2op(i,:)));
            curAdjIndx = tooCloseMis(wholeMatch,i,curAdjIndx);
            if ~isempty(curAdjIndx)
                pointAdjX{i}.adjPoi = unique([pointAdjX{i}.adjPoi curAdjIndx]);
                %             getPointSet(i,curAdjIndx,curYAdjIndex,wholeMatch);
                
                curMatch=wholeMatch;
                % multi-angle operation
                curMatch(curAdjIndx,4:5) = rotatePoi(curMatch(i,4:5),curMatch(curAdjIndx,4:5),rotationAngSet(j));
                locSimMea = calLocSim(i,curAdjIndx',curMatch,nameta1);
                locDisMea = calLocDis(i,curAdjIndx',curMatch);
                tmpPoiW(curAdjIndx) = tmpPoiW(curAdjIndx) + nameta*locSimMea+(1-nameta)*locDisMea;
                tmpPoiNum(curAdjIndx) = tmpPoiNum(curAdjIndx)+1; 
            end
        end
        adjIndx = find(tmpPoiW);
        w = tmpPoiW(adjIndx)./tmpPoiNum(adjIndx);
        pointAdjX{i}.simiVec = [pointAdjX{i}.simiVec adjIndx(find(w>theta))'];
    end
    pointAdjX{i}.simiVec = unique(pointAdjX{i}.simiVec);
            xl = pointAdjX{i}.simiVec;
            if length(pointAdjX{i}.simiVec)>=2
                curVecLabel = pointAdjX{i}.vecCluLabel;
                haveLabelIdx = find(curVecLabel);
                if isempty(haveLabelIdx)
                    if isempty(getAdjPoiLabel(pointAdjX,pointAdjX{i}.simiVec))
                        pointAdjX{i}.vecCluLabel = unique([pointAdjX{i}.vecCluLabel;labelNum]);
                        pointAdjX = setadjPoiLabel(pointAdjX,pointAdjX{i}.simiVec,labelNum);
                        labelNum=labelNum+1;
                    else
                        pointAdjX{i}.vecCluLabel=[pointAdjX{i}.vecCluLabel;getAdjPoiLabel(pointAdjX,pointAdjX{i}.simiVec)];
                    end
                else
                    %The current point has already been assigned a vector direction label, the similar matching vector is passed to it.
                    pointAdjX = setadjPoiLabel(pointAdjX,pointAdjX{i}.simiVec,curVecLabel);
                end
            end
        
    end

end

function labelSet = getAdjPoiLabel(pointInfo,index)
        labelSet = [];
        for i=1:length(index)
            labelSet = [labelSet;pointInfo{index(i)}.vecCluLabel];
        end
        labelSet=unique(labelSet);
end
%Assign the same label to similar neighboring points.
function pointInfo = setadjPoiLabel(pointInfo,index,labelNum)
for i=1:size(index)
    pointInfo{index(i)}.vecCluLabel = unique([pointInfo{index(i)}.vecCluLabel;labelNum]);
end
end
function findAdjIndx = tooCloseMis(wholeMatch,curIndx,curAdjIndx)
index = [curIndx curAdjIndx];
closeIndx=[];
for i=1:length(index)-1
    for j=i+1:length(index)
        if  sqrt((wholeMatch(index(i),1:2)-wholeMatch(index(j),1:2)).^2+(wholeMatch(index(i),4:5)-wholeMatch(index(j),4:5)).^2)<3
            closeIndx = [closeIndx;index(j)];
        end
    end
end
findAdjIndx = setdiff(curAdjIndx,closeIndx);
end

function pointSet = rotatePoi(srcPoi,srcAdjPoi,rotateAngle)

M1 = [1 0 srcPoi(1);
    0 1 srcPoi(2);
    0 0 1];
angles=rotateAngle;
A=[cosd(angles) -sind(angles) 0;
    sind(angles) cosd(angles) 0;
    0 0 1];
M2 = [1 0 -srcPoi(1);
    0 1 -srcPoi(2);
    0 0 1];

poix = M1*A*M2*[srcAdjPoi(:,1) srcAdjPoi(:,2) ones(size(srcAdjPoi,1),1)]';;
pointSet = poix(1:2,:)';
end