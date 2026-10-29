%DEMO 此处显示有关此函数的摘要
%   此处显示详细说明
clear;clc;
close all
%% set configuration
dataChosen.datasetLab='DAISY';
dataChosen.datasetName={'herzjesu' 'fountain' 'bearA' 'bearB' 'bearC' 'bearD' 'buildingA' 'buildingB'};
options.affineWarpFlag=false;
options.nameta=0.80;options.nameta1=0.80;options.theta=0.90;options.scale=1;
options.visible='on';
%%
im1=imread("./data/bearB0.jpg");
im2=imread("./data/bearB1.jpg");
if(size(im1,3))>1
    im1 = rgb2gray(im1);
end
if(size(im2,3))>1
    im2 = rgb2gray(im2);
end
% Some dataset like FIRE has its mask file for get ROI.
mask1=[];
mask2=[];
row1 = size(im1,2);
row2 = size(im2,2);

if ~isempty(mask1)
    im1 = getROI(im1,mask1);
    im2 = getROI(im2,mask2);
end
%% More initial match fuse and GroundTruth
load("./data/bearB01.mat");
wholeMatch =[X(:,[2 1]) (1:size(X,1))' Y(:,[2 1]) (1:size(Y,1))'];
trueMatch = [X(CorrectIndex,[2 1]) Y(CorrectIndex,[2 1])];
trueIdx=CorrectIndex;
inlierRate= size(trueMatch,1)*1.0./size(wholeMatch,1);
showmatch(im1,im2,wholeMatch(trueIdx,[1 2]),wholeMatch(trueIdx,[4 5]),0,options.visible,'Groundtruth matchs');
%%
[~,~,ib] = intersect(trueMatch,wholeMatch(:,[1 2 4 5]),'rows');
[matchIdx1] = findPoiMatchLabel(wholeMatch(:,[2 1 3 5 4 6]),options);
options.affineWarpFlag=true;
[matchIdx2] = findPoiMatchLabel(wholeMatch(:,[2 1 3 5 4 6]),options);
matchIdx=matchIdx1+matchIdx2;
%% get all method match
proMatch = [wholeMatch(find(matchIdx),[1 2 4 5]) find(matchIdx)];
showmatch(im1,im2,proMatch(:,1:2),proMatch(:,3:4).*options.scale,0,options.visible,'final matchs');