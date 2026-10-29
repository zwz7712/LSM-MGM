function [pointPairClusterLabelFin] = findPoiMatchLabel(match,option)
%FINDPOIMATCHLABEL6 此处显示有关此函数的摘要
%   此处显示详细说明
nameta=option.nameta;nameta1=option.nameta1;theta=option.theta;
matchNum = size(match,1);
rotationAngSet=[-70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70];
if option.affineWarpFlag
    wholeMatch=zeros(size(match));
    try
        wholeMatch(:,1:2) = affineL2e(match(:,1:2),match(:,4:5),1);
        wholeMatch(:,3:6)=match(:,3:6);
    catch
        pointPairClusterLabelFin=zeros(size(match,1),1);
        return
end
else
    wholeMatch=match;
end
%%
K=[10 15 20 25 30];K=20;
if matchNum<=6
    K=[1:1:matchNum-1];
else
    if matchNum<=30
        K=[mod(matchNum-1,5):5:matchNum-1];
    end
end
% pointSetX结构体：vecCluLabel:当前点的类别label; simiVec:领域内相似匹配向量索引；adjPoi:邻接点索引
[pointSetX,labelNum] = getPointDesV2(wholeMatch,K,nameta,rotationAngSet,nameta1,theta);
%% get match vector cluster label
operaIndx = (1:matchNum)';
sameCLuLabel=cell(labelNum,1);
% 获取每个点的聚类标签，如果多个标签取最小的一个，并传递给相同标签的其他点
for i=1:matchNum
    curLabel = pointSetX{i}.vecCluLabel;
    if length(curLabel)>1
        curMinLabel = min(curLabel);
        sameCLuLabel{curMinLabel}=setdiff(curLabel,curMinLabel);
    end
end
rootFlag=ones(size(sameCLuLabel,1),1);
pointPairClusterLabelFin=zeros(size(wholeMatch,1),1);
for i=1:labelNum
    if rootFlag(i)>0
        curClusTree = getFinClu(sameCLuLabel,i);
        rootFlag(i)=-i;
        rootFlag(curClusTree)=-i;
    end
end
for i=1:size(wholeMatch,1)
    if ~isempty(pointSetX{i}.vecCluLabel)
        pointPairClusterLabelFin(i)=-rootFlag(pointSetX{i}.vecCluLabel(1));
        operaIndx(i)=0;
    end
end

operaIndx = operaIndx(find(operaIndx));
for i=1:size(operaIndx,1)
    curPointSet = pointSetX{operaIndx(i)};
    simVecLabInd = find(pointPairClusterLabelFin(curPointSet.simiVec));
    % 存在相似结果属于不同类的情况
    if ~isempty(simVecLabInd)
        curSimVec = curPointSet.simiVec;
        labSet = unique(pointPairClusterLabelFin(curSimVec(simVecLabInd)),'rows');
        pointPairClusterLabelFin(operaIndx(i)) = labSet(1);
    end
end
% 消除中间不稳定的聚类结果
for i=1:labelNum
    curLabelPoiIdex = find(pointPairClusterLabelFin==i);
    if length(curLabelPoiIdex)<=2
        pointPairClusterLabelFin(curLabelPoiIdex)=0;
    end
    if size(unique(wholeMatch(curLabelPoiIdex,1:2),'rows'),1)<=2 || size(unique(wholeMatch(curLabelPoiIdex,4:5),'rows'),1)<=2
        pointPairClusterLabelFin(curLabelPoiIdex)=0;
    end
end
% 落到别的类簇位置的点进行判别，合并或消除更加不稳定的类？
unLabel = unique(pointPairClusterLabelFin);
for i=2:size(unLabel,1)
    curLabel = unLabel(i);
    curLbelPoiIndex = find(pointPairClusterLabelFin==curLabel);
    loc1 = wholeMatch(curLbelPoiIndex,1:2);
    loc2 = wholeMatch(curLbelPoiIndex,4:5);
    curLabelSize=size(loc1,1);
    try
        [ch1, v1] = convhull(loc1);
        [ch2, v2] = convhull(loc2);
        for j=1:size(pointPairClusterLabelFin,1)
            if pointPairClusterLabelFin(j)~=0 && pointPairClusterLabelFin(j)~=curLabel
                [~, v1t] = convhull([loc1;wholeMatch(j,1:2)]);
                [~, v2t] = convhull([loc2;wholeMatch(j,4:5)]);
                if v1t==v1 || v2t==v2
                    % 当前非该标签的点加入后凸包体积没变，即可认为点在该凸包内
                    curTmpLabel = pointPairClusterLabelFin(j);
                    curTmpLabelSize=length(find(pointPairClusterLabelFin==curTmpLabel));
                    adjIndex=0;
                    if v1t==v1
                        [~,adjIndex]=min(sqrt(sum((loc1-repmat(wholeMatch(j,1:2),curLabelSize,1)).^2,2)));
                        
                    else
                        [~,adjIndex]=min(sqrt(sum((loc2-repmat(wholeMatch(j,4:5),curLabelSize,1)).^2,2)));
                    end
                    %                 pointSetX{curLbelPoiIndex(adjIndex)}.simiVec
                    %                 pointSetX{j}.simiVec
                    if curTmpLabelSize<curLabelSize && size(pointSetX{j}.simiVec,1)<size(pointSetX{curLbelPoiIndex(adjIndex)}.simiVec,1)
                        pointPairClusterLabelFin([j;pointSetX{j}.simiVec])=0;
                    else if curTmpLabelSize>curLabelSize && size(pointSetX{j}.simiVec,1)>size(pointSetX{curLbelPoiIndex(adjIndex)}.simiVec,1)
                            pointPairClusterLabelFin(curLbelPoiIndex)=0;
                        else if curTmpLabelSize<curLabelSize  && size(pointSetX{j}.simiVec,1)>size(pointSetX{curLbelPoiIndex(adjIndex)}.simiVec,1)
                                ax=1;
                            else if curTmpLabelSize>curLabelSize  && size(pointSetX{j}.simiVec,1)<size(pointSetX{curLbelPoiIndex(adjIndex)}.simiVec,1)
                                    ax=2;
                                end
                            end
                        end
                    end
                end
            end
        end
    catch
        continue;
    end
end
end



    
    function cluResult = getFinClu(cLulabelSet,index)
        cluResult=[];
        curLabelSet = cLulabelSet{index};
        if ~isempty(curLabelSet)
            for i=1:length(curLabelSet)
                cluResult=[cluResult;curLabelSet(i)];
                cluResult=[cluResult;getFinClu(cLulabelSet,curLabelSet(i))];
            end
        end
    end


