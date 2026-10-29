function [localDes] = localTrendMeas(wholeMatch,srcIndex,srcAdjIndex,rotationAngSet,nameta,wThre,localDes,nbinsRa,nbinsTheta,thetaInd,radiusInd)
%LOCALTRENDMEAS 判断当前匹配与其周围匹配
%   此处显示详细说明
% figure(20);plot(wholeMatch(:,1),wholeMatch(:,2),'b+',wholeMatch(:,4),wholeMatch(:,5),'ro');
% nbinsTotal = nbinsRa*nbinsTheta;pointNum=size(wholeMatch,1);
% figure()
for j=1:size(rotationAngSet,2)
    curMatch = wholeMatch;
%     tmpDes = zeros(nbinsRa*nbinsTheta,1);
%     localDes{srcIndex,j} = zeros(pointNum,nbinsTotal);
    curMatch(srcAdjIndex,4:5) = rotatePoi(wholeMatch(srcIndex,4:5),wholeMatch(srcAdjIndex,4:5),rotationAngSet(j));
    
%     plot(wholeMatch(srcAdjIndex,4),curMatch(srcAdjIndex,5),'b+');
    locSimMea = calLocSim(srcIndex,srcAdjIndex',curMatch); 
    locDisMea = calLocDis(srcIndex,srcAdjIndex',curMatch);
    w = nameta*locSimMea+(1-nameta)*locDisMea;
%     tmpDes((thetaInd-1)*nbinsRa+radiusInd) = size(find(w>wThre),1);
%     if srcIndex==457
%         x=wholeMatch(srcIndex,1);y=wholeMatch(srcIndex,2);[r,~]=findMaxdistan(wholeMatch(srcIndex,1:2),wholeMatch(srcAdjIndex,1:2));
%         sx=wholeMatch([srcIndex srcAdjIndex],1:2);sy=wholeMatch([srcIndex srcAdjIndex],4:5);
%         sx2=curMatch([srcIndex srcAdjIndex],1:2);sy2=curMatch([srcIndex srcAdjIndex],4:5);
%         figure(21);plot(wholeMatch(:,1),wholeMatch(:,2),'b+',wholeMatch(:,4),wholeMatch(:,5),'ro');hold on;plot([sx(:,1) sy(:,1)]',[sx(:,2) sy(:,2)]','k-');rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1 1]);hold off;
%         figure(22);plot(curMatch(:,1),curMatch(:,2),'b+',curMatch(:,4),curMatch(:,5),'ro');hold on;plot([sx2(:,1) sy2(:,1)]',[sx2(:,2) sy2(:,2)]','k-');rectangle('Position',[x-r,y-r,2*r,2*r],'Curvature',[1 1]);hold off;
%     end
%     if size(find(w>wThre),1)~=0
%         al=1;
%     end
    localDes{j,1}(srcIndex,(thetaInd-1)*nbinsRa+radiusInd) = size(find(w>wThre),1); %对应j旋转角度的 对应区域描述子
    localDes{j,2}{srcIndex,:} = [localDes{j,2}{srcIndex,:} srcAdjIndex(w>wThre)];   % 对应j旋转角度的 相似匹配向量在点集A中索引
%     pause(10);
end                        
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

poix = M1*A*M2*[srcAdjPoi(:,1) srcAdjPoi(:,2) ones(size(srcAdjPoi,1),1)]';
% poix2 = M1*A*M2*[srcAdjPoi(1) srcAdjPoi(2) 1]';
pointSet = poix(1:2,:)';
end

function [maxdis,pointInd] = findMaxdistan(srcPoint,adjPoint)
tmpDis = sqrt(sum((repmat(srcPoint,size(adjPoint,1),1)-adjPoint).^2,2));
[maxdis,pointInd] = max(tmpDis);
end