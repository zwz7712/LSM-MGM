function [pointSetLabelFinal] = getPointLabelNew(point,totleSize,match,mean_dist_global,setFlag)
%GETGRIDLAY 获得点集中各点的标识，欧式空间为基础，在点集的最大最小范围内对点集进行均分（？是否能进行更合理的空间划分）
% 更多的参考的是点对匹配关系localmeasure还是点间距离？
%   此处显示详细说明
%% Kmeans聚类
% pointSetNum = size(pointSet,1);
% clusterNum = floor(pointSetNum/10);
% [idx,C] = kmeans(pointSet,clusterNum);
% pointSetLabel = idx;
%% 网格分类
% pointSetRangeX = [max(PointSet);min(PointSet)];
%% 极坐标bin分类
% global option
% 在极坐标下分成12个方向，最少分成5段距离进行参考，若其中前三段都没有合适结果当前方向放弃考虑，若断断续续包含结果和全部包含结果一起考虑进，加入下一段距离进行分析
% 按照一定顺序重排序，保证选择的点集与点匹配的对应关系
wThre = 0.95;
wholeMatch = match;
if setFlag ==1
    [pointSet,~,ic] = intersect(point,wholeMatch(:,1:2),'rows');
else
    [pointSet,~,ic] = intersect(point,wholeMatch(:,4:5),'rows');
end
wholeMatch = wholeMatch(ic,:);
r_inner=0.125;
nbins_theta=12;
r_outer=2;nbins_r=5;
pointNum=size(pointSet,1);
in_vec=1:pointNum;
% pointSet=2;
r_array=real(sqrt(dist2(pointSet,pointSet))); % real is needed to
                                          % prevent bug in Unix version
                                          % 取实部防止在Unix系统下出现bug
theta_array=atan2(pointSet(:,2)*ones(1,pointNum)-ones(pointNum,1)*pointSet(:,2)',pointSet(:,1)*ones(1,pointNum)-ones(pointNum,1)*pointSet(:,1)')';

% create joint (r,theta) histogram by binning r_array and
% theta_array
mean_dist_global=[];
% normalize distance by mean, ignoring outliers 使用mean值来归一化极半径
if isempty(mean_dist_global)
   tmp=r_array(in_vec,:);
   tmp=tmp(:,in_vec);
   mean_dist=mean(tmp(:));
end
r_array_n=r_array/mean_dist;

% use a log. scale for binning the distances
r_bin_edges=logspace(log10(r_inner),log10(r_outer),5);  % 等间距分的bin边缘
r_array_q=zeros(pointNum,pointNum);
for m=1:nbins_r
   r_array_q=r_array_q+(r_array_n<r_bin_edges(m));   
end
% fz=r_array_q>0; % flag all points inside outer boundary

% put all angles in [0,2pi) range
theta_array_2 = rem(rem(theta_array,2*pi)+2*pi,2*pi);
% quantize to a fixed set of angles (bin edges lie on 0,(2*pi)/k,...2*pi
theta_array_q = 1+floor(theta_array_2/(2*pi/nbins_theta));

nbins=nbins_theta*nbins_r;  % 所有块总量
localMeasure = cell(pointNum,1);
localMeasureArray = zeros(pointNum,pointNum);
% pointSetLabel = zeros(pointNum,1);
% 计算局部60个极坐标相关结果
for i=1:pointNum
    %     fzn = fz(i,:);
    curPoiFitMat = zeros(nbins_r,nbins_theta);
    for thInd = 1:nbins_theta
        thI = find(theta_array_q(i,:)==thInd);
        for arrInd = 1:nbins_r
            arrI = find(r_array_q(i,:)==6-arrInd);
            adjPoiInd = intersect(thI,arrI);
            adjPoiInd(adjPoiInd==i)=[];
            maintrendFlag = 0;
            maintrendInd = [];
            % 在当前bin中该点的对应点对是否占主流匹配趋势
            if adjPoiInd~=0
                % 局部的bin是否不适合于当前条件，局部可能有较多噪声
%                 [maintrendFlag,maintrendInd] = localTrendMeas([wholeMatch(:,1:3) zeros(pointNum,1) wholeMatch(:,4:6) zeros(pointNum,1)],i,adjPoiInd',0);
                calLocSim(i,adjPoiInd',[wholeMatch(:,1:3) zeros(pointNum,1) wholeMatch(:,4:6) zeros(pointNum,1)]);
                % 另一个思路，在这里尝试匹配？
            end
            if maintrendFlag
                curPoiFitMat(arrInd,thInd) = 1;
                localMeasureArray(i, maintrendInd) = 1;
                localMeasure{i} = curPoiFitMat;
                figure(1);
                %                  plot(wholeMatch(:,1),wholeMatch(:,2),'b+',wholeMatch(:,4),wholeMatch(:,5),'ro');
                plot([wholeMatch(:,1) wholeMatch(:,4)]',[wholeMatch(:,2) wholeMatch(:,5)]','r-');
                 hold on
                 
                 curMatch=[];
                 curMatch = [wholeMatch(i,:);wholeMatch(maintrendInd,:)];
                 plot(curMatch(1,1),curMatch(1,2),'b+',curMatch(1,4),curMatch(1,5),'ro');
                 plot([curMatch(:,1) curMatch(:,4)]',[curMatch(:,2) curMatch(:,5)]','k-');
                 hold off
            end
        end
    end
    
end
labelNum = 1;
pointSetLabel = zeros(pointNum,1);
relatLabel = cell(pointNum,1);
% labeltrendArray= [];
% 给所有点赋予多个label标签
% 部分点对的点在某些点的一个bin内 但该点的匹配与那个点匹配并不一致，仅作为近似的点块中内点存在
for i=1:size(localMeasure,1)
    if pointSetLabel(i)==0
        if ~isempty(localMeasure{i})
            [arrInd,theInd] = find(localMeasure{i}); % 确定其极坐标
            oldLabInd=[];addLabInd=[];minDist=0;nearestInd=0;
            addLabIndtmp=[];oldLabIndTmp=[];
            for indLen = 1:size(arrInd,1)
                % 以每行为一个点与当前所有点的关联结果，符合对应极坐标的对应点索引
                adjInd = intersect(find(r_array_q(i,:)==(6-arrInd(indLen))),find(theta_array_q(i,:)==theInd(indLen)));
                adjInd(adjInd==i)=[];
%                 find(localMeasureArray(i,:)==1)
%                 adjInd = intersect(adjInd,find(localMeasureArray(i,:)==1));
                if ~isempty(adjInd)
                    addFitInd = find(pointSetLabel(adjInd)==0);
                    
                    % 找周围计算的点中没有进行分label的点索引
                    if ~isempty(addFitInd)
                        addLabIndtmp = adjInd(addFitInd);
                        addLabInd=[addLabInd addLabIndtmp];
                    end
                    % 找周围计算的点中已经进行分label的点索引
                    oldLabIndTmp = setdiff(adjInd,addLabIndtmp);
                    oldLabInd = [oldLabInd oldLabIndTmp];
%                     if ~isempty(oldLabIndTmp)
%                         for tmpInd = 1:size(oldLabIndTmp,2)
%                             if r_array_q(i,oldLabIndTmp(tmpInd))>minDist
%                                 nearestInd=oldLabIndTmp(tmpInd);
%                                 minDist = r_array_q(i,oldLabIndTmp(tmpInd));
%                             end
%                         end
%                     end
                end
            end
%             labeltrendArray = [];
            if ~isempty(oldLabInd)
                oldLabCla = unique(pointSetLabel(oldLabInd));
                for oldLabClaInd = 1:size(oldLabCla,1)
                    oldCurLabInd = find(pointSetLabel==oldLabCla(oldLabClaInd));
                    w = calLocSim(i,oldCurLabInd,[wholeMatch(:,1:3) zeros(pointNum,1) wholeMatch(:,4:6) zeros(pointNum,1)]);
                    if length(find(w>wThre))-1>=(length(oldCurLabInd)-1)*2/3.0
                        pointSetLabel(i) = oldLabCla(oldLabClaInd);
                        relatLabel{i} = [relatLabel{i};oldLabCla(oldLabClaInd)];
                        relatLabel{i} = unique(relatLabel{i});
                    end
                end
                if pointSetLabel(i)==0
                    pointSetLabel(i) = labelNum;
                    relatLabel{i} = [relatLabel{i};labelNum];
                    relatLabel{i}=unique(relatLabel{i});
                    labelNum=labelNum+1;
                end
            else
%                 curLabel=labelNum;
                pointSetLabel(i) = labelNum;
                relatLabel{i} = [relatLabel{i};labelNum];
                relatLabel{i}=unique(relatLabel{i});
%                 labeltrendArray{labelNum} = [labeltrendArray{labelNum}; (find(pointSetLabel==labelNum))']
                labelNum=labelNum+1;
            end
            
            if ~isempty(addLabInd)
                for addlabIndlen =  1:size(addLabInd,2)
                    addLabClaInd = find(pointSetLabel==pointSetLabel(i));
                    w = calLocSim(addLabInd(addlabIndlen),addLabClaInd,[wholeMatch(:,1:3) zeros(pointNum,1) wholeMatch(:,4:6) zeros(pointNum,1)]);
                    if length(find(w>wThre))-1>=(length(addLabClaInd)-1)*2/3.0
                        pointSetLabel(addLabInd(addlabIndlen))=pointSetLabel(i);
                        relatLabel{addLabInd(addlabIndlen)}=[relatLabel{addLabInd(addlabIndlen)};pointSetLabel(i)];
                        relatLabel{addLabInd(addlabIndlen)} = unique(relatLabel{addLabInd(addlabIndlen)});
%                         if i==32
%                             showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,4:5),pointSetLabel,1,2,['111']);
%                         end
                    end
                end
            end
            
            showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,4:5),pointSetLabel,1,2,['111']);
            a12=1;
        end
    else
        if ~isempty(localMeasure{i})
            [arrInd,theInd] = find(localMeasure{i});
            oldLabInd=[];addLabInd=[];minDist=0;nearestInd=0;
            addLabIndtmp=[];oldLabIndTmp=[];
            for indLen = 1:size(arrInd,1)
                adjInd = intersect(find(r_array_q(i,:)==(6-arrInd(indLen))),find(theta_array_q(i,:)==theInd(indLen)));
                adjInd(adjInd==i)=[];
%                 find(localMeasureArray(i,:)==1)
%                 adjInd = intersect(adjInd,find(localMeasureArray(i,:)==1));
                if ~isempty(adjInd)
                    addFitInd = find(pointSetLabel(adjInd)==0);
                    % 找周围计算的点中没有进行分label的点索引
                    if ~isempty(addFitInd)
                        addLabIndtmp = adjInd(addFitInd);
                        addLabInd=[addLabInd addLabIndtmp];
                    end
                    % 找周围计算的点中已经进行分label的点索引
                    oldLabIndTmp = setdiff(adjInd,addLabIndtmp);
                    oldLabInd = [oldLabInd oldLabIndTmp];
                    if ~isempty(oldLabIndTmp)
                        for tmpInd = 1:size(oldLabIndTmp,2)
                            if r_array_q(i,oldLabIndTmp(tmpInd))>minDist
                                nearestInd=oldLabIndTmp(tmpInd);
                                minDist = r_array_q(i,oldLabIndTmp(tmpInd));
                            end
                        end
                    end
                end
            end
            if ~isempty(oldLabInd)
                 oldLabCla = unique(pointSetLabel(oldLabInd));
                 for oldLabClaInd = 1:size(oldLabCla,1) 
                     oldCurLabInd = find(pointSetLabel==oldLabCla(oldLabClaInd));
                     w = calLocSim(i,oldCurLabInd,[wholeMatch(:,1:3) zeros(pointNum,1) wholeMatch(:,4:6) zeros(pointNum,1)]);
                     if length(find(w>wThre))-1>=(length(oldCurLabInd)-1)*2/3.0
                         pointSetLabel(i) = oldLabCla(oldLabClaInd);
                         relatLabel{i} = [relatLabel{i};oldLabCla(oldLabClaInd)];
                         relatLabel{i} = unique(relatLabel{i});
                         %                     pointSetLabel(i) = pointSetLabel(oldLabInd(1));
                     end
                 end
            end
            if ~isempty(addLabInd)
                for addlabIndlen =  1:size(addLabInd,2)
                    addLabClaInd = find(pointSetLabel==pointSetLabel(i));
                    w = calLocSim(addLabInd(addlabIndlen),addLabClaInd,[wholeMatch(:,1:3) zeros(pointNum,1) wholeMatch(:,4:6) zeros(pointNum,1)]);
                    if length(find(w>wThre))-1>=(length(addLabClaInd)-1)*2/3.0
                        pointSetLabel(addLabInd(addlabIndlen))=pointSetLabel(i);
                        relatLabel{addLabInd(addlabIndlen)}=[relatLabel{addLabInd(addlabIndlen)};pointSetLabel(i)];
                        relatLabel{addLabInd(addlabIndlen)} = unique(relatLabel{addLabInd(addlabIndlen)});
                    end
                end
            end
             showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,4:5),pointSetLabel,1,2,['111']);
             a12=1;
        end
            
%         end
       
    end
%     tmpPointSetLabel = zeros(pointNum,1);
%     tmpPointSetLabel(i)=1;
    
end
% 集合标签，得到最终分类结果(有问题，有着多个标签的结果不能直接集中到一个类别中，可能是中间的某些点符合多个类别标签，但几个类别区别较大)
labFlag = zeros(size(pointSetLabel,1),1);

labelNum = 1;
% 整合相同的点集标签
pointSetLabelFinal = getLabelInte(pointSetLabel,relatLabel);
showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,4:5),pointSetLabelFinal,1,3,['111']);
adjind=[];
for i=1:size(localMeasure,1)
    if ~isempty(localMeasure{i})
        adjind(i) = 1;
    else
        adjind(i) = 0;
    end
end
adjind=logical(adjind);
pointSetLabelFinal(~adjind,:)=0;
% option = figIndSet(option,'Mv2FxMatch');
showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,4:5),pointSetLabelFinal,1,4,['111']);
showcluster(pointSet,pointSetLabelFinal,4);
al = setdiff([find(pointSetLabelFinal==0);find(pointSetLabel==0)],find(pointSetLabel==0));
maxLabel = max(pointSetLabelFinal);
% noLabelInd = find(pointSetLabelFinal==0);
% pointSetLabelFinal(noLabelInd)= maxLabel+1:maxLabel+noLabelInd;
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
                curAdjlabel = [curAdjlabel;relatLabel{curInd(j)}];
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

function showcluster(pointset, label,figureIdx)
colorLib = ['r','g', 'b' 'c','m','y','k'];
labelInd = unique(label((find(label>0))));
len = length(labelInd);
figure(figureIdx);
hold on;
for i=1:size(labelInd,1)
    indx = find(label==labelInd(i));
    idx = mod(i,7);
    if idx == 0
        idx = 5;
    end
    plot(pointset(indx,1),pointset(indx,2),[colorLib(idx) 'o']); 
end
hold off
end