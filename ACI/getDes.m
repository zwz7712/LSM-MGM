function [localMeasure] = getDes(pointSet,wholeMatch,mean_dist_global)
%GETDES 此处显示有关此函数的摘要
%   此处显示详细说明
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
% mean_dist_global=[];
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
            % 在当前bin中该点的对应点对是否占主流匹配趋势
            if adjPoiInd~=0
                [maintrendFlag,~] = localTrendMeas(wholeMatch,i,adjPoiInd',0);
            end
            if maintrendFlag
                curPoiFitMat(arrInd,thInd) = 1;
                localMeasure{i} = curPoiFitMat;
            end
        end
    end
    
end
end

