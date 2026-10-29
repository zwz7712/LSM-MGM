function localDis = calLocDis(curIndex,curAdjIndex,wholeMatch)
%CAL 此处显示有关此函数的摘要
%   此处显示详细说明
curNum = size(curAdjIndex,1);
curPoint = wholeMatch(curIndex,1:2);
curCorPoint = wholeMatch(curIndex,4:5);
curAdjPoint = wholeMatch(curAdjIndex,1:2);
curAdjCorPoint = wholeMatch(curAdjIndex,4:5);
thet=15;
al = ((curAdjPoint-repmat(curPoint,curNum,1))-(curAdjCorPoint-repmat(curCorPoint,curNum,1))).^2;
x1=sqrt(sum((curAdjPoint-repmat(curPoint,curNum,1)).^2,2));
x2=sqrt(sum(((curAdjCorPoint-repmat(curCorPoint,curNum,1))).^2,2));
localDis = exp(-(sum((sqrt(sum((curAdjPoint-repmat(curPoint,curNum,1)).^2,2))-sqrt(sum(((curAdjCorPoint-repmat(curCorPoint,curNum,1))).^2,2))).^2,2))./(2*thet^2));
end

