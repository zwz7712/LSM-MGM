function [localSimilar] = calLocSim(curIndex,curAdjIndex,wholeMatch,l1)
%CALLOCSIM 根据提供的当前点对和其邻近点对计算局部相似性值
%   curIndex: 当前匹配在整体匹配集合中的索引
%   curAdjIndex: 当前匹配点对的邻域点对索引（在点集A中）

% l1 = 0.5; % 长度惩罚系数？
l2 = 1-l1; % 角度惩罚系数？
thet1 = 0.4; % 距离高斯惩罚系数
thet2 = 0.8; % 角度高斯惩罚系数
curNum = size(curAdjIndex,1);
curPoint = wholeMatch(curIndex,1:2);
curCorPoint = wholeMatch(curIndex,4:5);
% adjPoi = wholeMatch(curAdjIndex,:);
curVec = curCorPoint - curPoint; % current translation vector
curtrendSet = wholeMatch(curAdjIndex,4:5) - wholeMatch(curAdjIndex,1:2);
curtrendNum = size(curtrendSet,1);
simIdx = find(sqrt(repmat(sum(curVec.^2,2),curtrendNum,1)+sum(curtrendSet.^2,2))<=5);
localSimilar = zeros(curtrendNum,1);
if ~isempty(simIdx)
    curtrendSimSet = curtrendSet(simIdx,:);
    tmplocalSimilar1 = calLocDis(curIndex,curAdjIndex(simIdx),wholeMatch);
    localSimilar(simIdx) = tmplocalSimilar1;
    tmplocalSimilar2 = [];
    if size(simIdx,2)~=curtrendNum
        diffSimIdx = setdiff(1:curtrendNum,simIdx);
        tmplocalSimilar2 = getSimVal(curVec,curtrendSet(diffSimIdx,:),l1,l2,thet1,thet2,size(diffSimIdx,2));
        localSimilar(diffSimIdx) = tmplocalSimilar2;
    end
    
else
    localSimilar=getSimVal(curVec,curtrendSet,l1,l2,thet1,thet2,curNum);
end
end

function vec = adjustVec(vec)
idx = find(sum(vec.^2,2)<=1);
vec(idx,:)=zeros(size(idx,1),2);
end

function localSimilar = getSimVal(curVec,curtrendSet,l1,l2,thet1,thet2,curNum)
% fitgeotrans
angles = repmat(atan2d(curVec(2),curVec(1)),size(curtrendSet,1),1)-atan2d(curtrendSet(:,2),curtrendSet(:,1));
ind = find(angles>180);
angles(ind)=angles(ind)-360*sign(angles(ind));

a=repmat(sqrt(sum(curVec.^2,2)),curNum,1); 
b=sqrt(sum(curtrendSet.^2,2)); 

localSimilar=zeros(curNum,1);
if a==0
    localSimilar(find(b==0))=1;
    localSimilar(find(b))=0;
else
    disPnis = max([a,b],[],2)./min([a,b],[],2); 
    sv1 = l1*exp(-(disPnis-1).^2./(2*thet1^2));
    sv2 = l2*exp(-abs(angles.*disPnis/180*pi).^2./(2*thet2^2));
    localSimilar =sv1+sv2;
end
end