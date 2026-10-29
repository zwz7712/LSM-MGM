function [ pointIndex,adjIndex] = adjacentPoint(pointSet,cell,index )
% 功能：在wholeMatch中寻找以该cell为中心的邻接的点
% Input:
% pointSet 1:2维是点二维坐标，3:4维是对应cell索引
% cell 当前选择的cell
% index 当前点在set中的索引位置
% Output:
% pointIndex 选中的点的二维坐标
% adjIndex 选中点在点集中的索引位置
pointIndex = [];
adjIndex = [];
set = pointSet;
set(index,3) = -1;
num = size(set,1);

% 搜索包含本cell在内的点索引
curCell = cell;
indexS = find(sum(abs(set(:,3)-repmat(curCell,num,1)),2)==0);
if(~isempty(indexS))
    pointIndex = set(indexS,1:2);
    adjIndex = indexS;
end

end