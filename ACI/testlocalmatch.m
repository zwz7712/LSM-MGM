function [outputArg1,outputArg2] = testlocalmatch(movPoi,fixPoi,wholeMatch)
%TESTLOCALMATCH 此处显示有关此函数的摘要
%   此处显示详细说明
%% 极坐标bin分类
% 在极坐标下分成12个方向，最少分成5段距离进行参考，若其中前三段都没有合适结果当前方向放弃考虑，若断断续续包含结果和全部包含结果一起考虑进，加入下一段距离进行分析
r_array=real(sqrt(dist2(fixPoi(:,1:2),fixPoi(:,1:2))));

tmp=r_array(in_vec,:);
tmp=tmp(:,in_vec);
mean_dist=mean(tmp(:));
mean_dist_array = [mean_dist mean_dist*1.5 mean_dist*2 mean_dist*3 mean_dist*2/3 mean_dist/2 mean_dist/3];
   
des1=getDes(movPoi,wholeMatch,[]);
des2Cur=[];
match=[];
for i=1:size(mean_dist_array)
des2Cur =  getDes(fixPoi,wholeMatch,mean_dist_array(i));
match = 
end
end

