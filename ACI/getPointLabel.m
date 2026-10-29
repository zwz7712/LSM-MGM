function [pointPairClusterLabelFinal] = getPointLabel(point,totleSize,match,mean_dist_global,setFlag)
%GETPOINTLABEL 此处显示有关此函数的摘要
%   此处显示详细说明
wThre = 0.95;
wholeMatch = match;
pointSet=[];
if setFlag ==1
    [pointSet,ip,iw] = intersect(point,wholeMatch(:,1:2),'rows');
else
    [pointSet,ip,iw] = intersect(point,wholeMatch(:,4:5),'rows');
end
pointNum = size(pointSet,1);
pointSetLabel=zeros(length(point),1);
notFitLabel=cell(length(point),1);
wholeMatch = wholeMatch(iw,:);
PoiRelaLab = cell(pointNum,1);
trendSet=[];
labelNum = 1;
% 1、只考虑局部一致性的分类
for pointInd = 1:pointNum
    fitInd=[];
    curAdjInd = [1:pointNum];
    curAdjInd(find(curAdjInd==pointInd))=[];
    w = calLocSim(pointInd,curAdjInd',[wholeMatch(:,1:3) zeros(pointNum,1) wholeMatch(:,4:6) zeros(pointNum,1)]);
    fitInd = find(w>=wThre);
    curAdjLabSet=[];
    for i=1:size(fitInd,2)
        curAdjLabSet = [curAdjLabSet PoiRelaLab{curAdjInd(fitInd)}];
    end
    curAdjLabSet=unique(curAdjLabSet);
    if isempty(find(curAdjLabSet))
        %        curLabelMean = mean(w(fitInd));
        %        pointSetLabel([pointInd;curAdjInd(fitInd)])=labelNum;
        PoiRelaLab = valCell(PoiRelaLab,[pointInd curAdjInd(fitInd)],labelNum);
        
        notFitLabel{labelNum} = [notFitLabel{labelNum} curAdjInd(find(w<wThre))];notFitLabel{labelNum}=unique(notFitLabel{labelNum});
        trendSet{labelNum} = [pointInd curAdjInd(fitInd)];trendSet{labelNum}=unique(trendSet{labelNum});
        labelNum=labelNum+1;
    else
%         labelCal = unique(pointSetLabel(fitInd));
        for labCalInd = 1:size(curAdjLabSet,2)
            curLabel = curAdjLabSet(labCalInd);
            w = calLocSim(pointInd,trendSet{curLabel}',[wholeMatch(:,1:3) zeros(pointNum,1) wholeMatch(:,4:6) zeros(pointNum,1)]);
            if length(find(w>=wThre))>=length(w)*2/3
                trendSet{curLabel} = [trendSet{curLabel} pointInd];trendSet{curLabel}=unique(trendSet{curLabel});
                PoiRelaLab{pointInd} = [PoiRelaLab{pointInd} curLabel];PoiRelaLab{pointInd} =unique(PoiRelaLab{pointInd});
                % 针对当前在判断的点，对于当前比较的点块中不太相似的点该如何处理
%                 notFitLabel{curLabel} = [notFitLabel{labelNum};curAdjInd(find(w<wThre))];
            else
                 notFitLabel{curLabel} = [notFitLabel{curLabel} pointInd];notFitLabel{curLabel}=unique(notFitLabel{curLabel});
            end
        end
        % 相关的旧类别都和本点的匹配趋势不一致
        if isempty(PoiRelaLab{pointInd})
            PoiRelaLab = valCell(PoiRelaLab,[pointInd curAdjInd(fitInd)],labelNum);
            notFitLabel{labelNum} = [notFitLabel{labelNum} curAdjInd(find(w<wThre))];notFitLabel{labelNum}=unique(notFitLabel{labelNum});
            trendSet{labelNum} = [pointInd curAdjInd(fitInd)];trendSet{labelNum}=unique(trendSet{labelNum});
            labelNum=labelNum+1;
        end
    end
end
% 考虑了距离的重分类
maxLen=1;
for i=1:pointNum
    len = length(PoiRelaLab{i});
    if len>maxLen
    maxLen = len;
    end
end
pointPairClusterLabel=zeros(pointNum,1);
for i=1:1
    
    for j=1:pointNum
        curlabset = PoiRelaLab{j};
         len = length(curlabset);
    if len>=i
        pointPairClusterLabel(j) = curlabset(i);
    end
    end
%     showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,4:5),pointPairClusterLabel,1,i,'111')
end

pointPairClusterLabel = getLabelInte(pointPairClusterLabel,PoiRelaLab);
% showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,4:5),pointPairClusterLabelFinal,1,i,'111')
pointPairClusterLabelFinal = zeros(size(point,1),1);
% [~,~,icv] = intersect( point,pointset);
for i=1:size(point,1)
    ind=0;
    ind = find(sum(abs(pointSet-repmat(point(i,:),size(pointSet,1),1)),2)==0);
    if ind
        pointPairClusterLabelFinal(i) = pointPairClusterLabel(ind);
    end
end
% pointPairClusterLabelFinal = pointPairClusterLabel();
end

function cell = valCell(cellset,ind,val)
cell=cellset;
for  i = 1:length(ind)
    cell{ind(i)} = [cell{ind(i)} val];
    cell{ind(i)}=unique(cell{ind(i)});
end
end

function [pointSetLabelFinal] = getLabelInte(pointSetLabel,relatLabel)
pointNum = size(pointSetLabel,1);
pointSetLabelFinal = zeros(size(pointSetLabel,1),1);
labNum=1;
% intersect(find(pointSetLabelFinal));
labelArr = unique(pointSetLabel);
for i=1:pointNum
    if pointSetLabel(i)~=0
        if pointSetLabelFinal(i)==0
            curLab = pointSetLabel(i);
            curInd = find(pointSetLabel==curLab);
            curLabNum = labNum;
            %         pointSetLabelFinal(curInd)=labNum;
            curAdjlabel=[];
            for j=1:size(curInd,1)
                curAdjlabel = [curAdjlabel;relatLabel{curInd(j)}'];
            end
            curAdjlabel = setdiff(unique(curAdjlabel),curLab);
            arrayLabel = [];
            arrayInd = [];
            noLabeArrInd = [];
            if ~isempty(curAdjlabel)
                for j=1:size(curAdjlabel,1)
                    ind = find(pointSetLabel==curAdjlabel(j));
                    if pointSetLabelFinal(ind)~=0
                        arrayLabel = [arrayLabel;pointSetLabelFinal(ind)];
                        arrayInd = [arrayInd;ind];
                    else
                        noLabeArrInd = [noLabeArrInd;ind];
                    end
                end
            end
            if ~isempty(arrayLabel)
                curLabNum = min(arrayLabel);
                % 根据本点集相关label的点修改编号
                
                %             arrayInd
            else
                
                labNum=labNum+1;
            end
%             if ~isempty(find(curInd==66)) || ~isempty(find(arrayInd==66)) || ~isempty(find(noLabeArrInd==66))
%                 a=1
%             end
            % 当前及同样label的编号增加
            pointSetLabelFinal(curInd) = curLabNum;
            pointSetLabelFinal(i) = curLabNum;
            % 相关联label的编号增加或修改
            pointSetLabelFinal(arrayInd) = curLabNum;
            pointSetLabelFinal(noLabeArrInd) = curLabNum;
        end
    end
end
end