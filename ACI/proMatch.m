   function [ finalMatch,finalPointPairClusterLabel ] = proMatch(wholeMatch,totalSize,im1,im2,trueMatch,option)
%初匹配的后处理
% 配准后的结果点对，其相邻的点对其对应点也应该在另一幅图的相邻区域
% adjacent cell invariance(ACI)限制 (选取到外点后迭代处理)
% 利用该限制，对前述匹配结果进行后处理，并构建最初的点块模型
%Input:
%   wholematch:输入初匹配点对，N x D N个点对， D维：1:2 点1xy坐标，3 所在cluster索引 4 原点集编号， 5:6 点2xy坐标 7 所在cluster索引 8 原点集编号
%   cellSize 一个cell的尺寸 1 x 2 row col
%   totalSize 数据点集的长宽或图片大小 2 x 2 row col
%Output:
%   finalMatch:     最终正确匹配对
%   finalPointPairClusterLabel: 初始匹配对对应cluster分类
%
global options
a = sqrt(sum((wholeMatch(:,1:2)-wholeMatch(:,4:5)).^2,2));
num = size(wholeMatch,1);
pointPairClusterLabel = zeros(num,1); % 匹配点对聚类标签
nameta=0.2;
%% 以点集1为基础点集，根据当前点对在点集1的点于其他点对在点集1的点是否与其邻近做搜寻标准，找出这些点对在点集2的点是否也是和该点对在点集2的点存在一致
% 过于强调局部的一致性，对于在全局上考虑上不合适的点（在局部不存在较多正确匹配点对，而碰巧存在较多相似的错误结果，需要平滑函数来限制外点防止过拟合）
% [pointPairClusterLabel] = findPoiMatchLabel2(pointPairClusterLabel,wholeMatch,totalSize,1:num);
[~,~,ib] = intersect(trueMatch,wholeMatch(:,[1 2 4 5]),'rows');
showMatch(wholeMatch(:,[2 1]),wholeMatch(:,[5 4]),ib,ib,7,['GT Num:' num2str(length(ib))],1);
[match] = findPoiMatchLabel4(wholeMatch(:,[2 1 3 5 4 6]),1:num,ib,nameta);
% [match] = findPoiMatchLabel5(wholeMatch(:,[2 1 3 5 4 6]),1:num,ib,nameta);
showMatch(wholeMatch(:,[2 1]),wholeMatch(:,[5 4]),find(match),find(match),6,['test Num:' num2str(length(find(match)))]);
finalPointPairClusterLabel = match;

pointPairClusterLabel = zeros(num,6);
% getClus(wholeMatch,find(match));
% 根据去误匹配结果聚类
haveMatch = find(match);
tmpLabel = DBSCANClus(wholeMatch(:,[2 1 3 5 4 6]),find(match~=0),nameta);
% 根据初匹配结果聚类
fixedIdx = find(match~=0 & tmpLabel~=0);
pointPairClusterLabel(fixedIdx) = tmpLabel(fixedIdx);

% a = setdiff(ib,find(pointPairClusterLabel),'rows');
showMatch(wholeMatch(:,[2 1]),wholeMatch(:,[5 4]),find(match),find(match),6,['test Num:' num2str(length(find(match)))]);
% showMatch(wholeMatch(:,[2 1]),wholeMatch(:,[5 4]),ib,ib,7,['GT Num:' num2str(length(ib))]);
% showMatchLabel(wholeMatch(:,[2 1]),wholeMatch(:,[5 4]),tmpLabel,option.visible,8,'testLabel');
len = length(unique(pointPairClusterLabel((find(pointPairClusterLabel>0)))));
matchNum = length(find(pointPairClusterLabel>0));
if(len==0)
    error('find no right match in this pair');
else
    option = figIndSet(option,'Mv2FxMatch');
    showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,5:6),pointPairClusterLabel,option.visible,option.figureInd.Mv2FxMatch,['The first Match from moving to fixing. MatchNum: ' num2str(matchNum)]);
    %     showMatchLabel(wholeMatch,pointPairClusterLabel,display_flag,2,'The first Match from moving to fixing');
end
% A=1;
%% 反向检索
pointPairClusterLabelTmp = zeros(num,1);
pointPairClusterLabelInv = zeros(num,1);
% 反向检索的完整
tempwholeMatch = zeros(size(wholeMatch,1),size(wholeMatch,2));
tempwholeMatch(:,5:8) = wholeMatch(:,1:4);
tempwholeMatch(:,1:4) = wholeMatch(:,5:8);
tempTotalSize = zeros(size(totalSize));
tempTotalSize(1,:) = totalSize(2,:);
tempTotalSize(2,:) = totalSize(1,:);
pointPairClusterLabelInv = findPoiMatchLabel2(pointPairClusterLabelInv,tempwholeMatch,tempTotalSize,1:num);
len = length(unique(pointPairClusterLabelInv((find(pointPairClusterLabelInv>0)))));
matchNum = length(find(pointPairClusterLabelInv>0));
if(len==0)
    error('find no right match in this pair');
else
    option = figIndSet(option,'Fx2MvMatch');
    showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,5:6),pointPairClusterLabelInv,option.visible,option.figureInd.Fx2MvMatch,['The Second Match from moving to fixing. MatchNum: ' num2str(matchNum)]);
    %     showMatchLabel(wholeMatch,pointPairClusterLabel,display_flag,2,'The first Match from moving to fixing');
end
%% Get stable matches %%%%% Bayesian Integration place it ???
% matchNum = length(find(pointPairClusterLabelInv));
% showMatchLabel(tempMatch(:,1:2),tempMatch(:,5:6),pointPairClusterLabelTmp,display_flag,3,['The Second Match from fixing to Matching. MatchNum: ' num2str(matchNum)]);
label = unique(pointPairClusterLabel(find(pointPairClusterLabel>0)));
labelLength = length(label);
labelCor = cell(labelLength,2); %dim 1: label index, match index; dim 2: label index in inverse match, match index
% curPoiLabel = cell(1,2);
for i=1:num
    % 主要参考moving set label？进行label redesign
    if(pointPairClusterLabel(i)>0)
        labelInd = find(label==pointPairClusterLabel(i));
        labelCor{labelInd,1} = [labelCor{labelInd,1};label(labelInd), i];
    end
    if(pointPairClusterLabel(i)>0 && pointPairClusterLabelInv(i)>0)
        % 结合2次匹配结果考虑正反向对应label的参考
        labelCor{labelInd,2} = [labelCor{labelInd,2};pointPairClusterLabelInv(i), i];
        %
        %         pointPairClusterLabelTmp(i) = pointPairClusterLabelInv(i);
        temp = zeros(num,1);tempInv = temp;
        temp(i) = pointPairClusterLabel(i);tempInv(i) = pointPairClusterLabelInv(i);
%     else
%         % 正反向匹配有一次没有匹配上，算作不稳定的匹配结果
    end
end

finLabNum=1;
for i=1:labelLength
    curMatch = labelCor{i,1};
    curMatchInv = labelCor{i,2};
    if ~isempty(curMatchInv)
        curInvLabel = unique(curMatchInv(:,1));
        sameLabelLen = length(curInvLabel);
        if sameLabelLen==1
            % 当前正方向的正确匹配结果的同一label 对应于反方向label也只有一种结果
            pointPairClusterLabelTmp(curMatch(:,2))=finLabNum;
            finLabNum=finLabNum+1;
        else 
            % 当前正方向的正确匹配结果的同一label 对应于反方向label有多种结果
            % 根据第一个匹配的label进行合并？
            for j=1:sameLabelLen
                tmpIdx = find(curMatchInv(:,1)==curInvLabel(j));
                pointPairClusterLabelTmp(curMatchInv(tmpIdx,2)) = finLabNum;
            end
        end
        finLabNum=finLabNum+1;
    else
        % 当前正方向的正确匹配结果的同一label 没有反方向label？？？
    end
end

option = figIndSet(option,'BBFMatch');
matchNum = length(find(pointPairClusterLabelTmp>0));
showMatchLabel(wholeMatch(:,1:2),wholeMatch(:,5:6),pointPairClusterLabelTmp,option.visible,option.figureInd.BBFMatch,['BBF Merge Match. MatchNum:' num2str(matchNum)]);

%% 汇总最终结果，对相同区域、相同移动趋势的不同簇结果进行合并，同一簇不同移动趋势的结果进行分割
% 对于第一个各个cluster的label，存在情况：
%       1.在第二个label中存在多个label对应于某一个label的情况;
%       2.一对一的情况（稳定配对，直接保留;基本为外点的情况？;多个属于第一个label对应于一个第二个label的情况？）
%       3.无对应情况 （不稳定配对，直接消除）

%衡量向量是否属于一个簇（向量长度，夹角，于图中的位置）


finalMatch = wholeMatch(find(pointPairClusterLabelTmp>0),:);
fitIdx = find(pointPairClusterLabelTmp>0);
% % 最小二乘？ RANSAC
inlierNum = size(finalMatch,1);
if inlierNum>=8
    [fRANSAC,inliers] =  estimateFundamentalMatrix(finalMatch(:,1:2),finalMatch(:,5:6),'Method','RANSAC','NumTrials',2000);
else
    x1 = finalMatch(:,1);
    y1 = finalMatch(:,2);
    x2 = finalMatch(:,5);
    y2 = finalMatch(:,6);
    % 匹配点对多于6少于8
    if inlierNum>=6
        try
            t_fundus = fitgeotrans([x2 y2],[x1 y1],'polynomial',2);
            transformedIm = transformPointsInverse(t_fundus,[x1 y1]);
            dist = sqrt(sum((transformedIm-[x2 y2]).^2,2));
            [Va,idx] = sort(dist);
            % 清洗后内点选择：1、数量前20%会不会删了太多了，这样保留的结果很通过点对簇标识难传播开
            %                2、最小距离的150%的比例内范围
            %                 inliers = idx(1:floor(corNum*0.2));
            indx = idx(1:find(Va<Va(1)*2));
            inliers(indx) = 1;
        catch ErrorInfo
            inliers = ones(inlierNum,1);
        end
    else
        % 匹配点对多于3少于6
        if inlierNum>=3
            try
                t_fundus = fitgeotrans([x2 y2],[x1 y1],'affine');
                transformedIm = transformPointsInverse(t_fundus,[x1 y1]);
                dist = sqrt(sum((transformedIm-[x2 y2]).^2,2));
                [Va,idx] = sort(dist);
                % 清洗后内点选择：1、数量前20%会不会删了太多了，这样保留的结果很通过点对簇标识难传播开
                %                2、最小距离的150%的比例内范围 （错误，如果最近点就在该点旁边，那后续所有可参考点都变得没有意义）
                %                     inliers = idx(1:floor(corNum*0.2));
                indx = idx(1:find(Va<Va(1)*1.5));
                inliers(indx) = 1;
            catch ErrorInfo
                inliers = ones(inlierNum,1);
            end
        else
            % 匹配点对少于3
            % 单匹配保留， 多匹配非同方向全部消除
            if length(ind)==1
                inliers = 1;
            else
                [~,isSameDire] = findMainTrend(ind(1),ind(2),wholeMatch(:,1:2),wholeMatch(:,5:6));
                if isSameDire
                    inliers = ones(inlierNum,1);
                else
                    inliers = [];
                end
            end
        end
    end
end
% if is emp
finalMatch = finalMatch(inliers,:);
finalPointPairClusterLabel = zeros(num,1);
finalPointPairClusterLabel(fitIdx(inliers)) = pointPairClusterLabelTmp(fitIdx(inliers));
% idxMatch = idxMatch(inliers,:);
end

function [angle] = motionAngle(matchPoi1,matchPoi2)
curVec = matchPoi2-matchPoi1;
angle = repmat(atan2d(curVec(1),curVec(2)),size(adjPoiTrend,1),1)-atan2d(adjPoiTrend(:,1),adjPoiTrend(:,2));
ind = find(angle>180);
if(angle>180)
    angle=angle-360*sign(angle);
end
end
