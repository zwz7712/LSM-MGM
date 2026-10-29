function showmatch(I1_p,I2_p,loc1,loc2,s,visible,titles)
% r 非对称高斯参数
% Input：
%       I1_p 输入图像1
%       I2_p 输入图像2
%       loc1 N x 2维 对应点对在X点集中的点坐标 第一维是总坐标 第二维横坐标
%       loc2 N x 2维 对应点对在Y点集中的点坐标

if(strcmp(visible,'on'))
    if(size(I1_p,3))>1
        I1 = I1_p(s+1:end-s,s+1:end-s,:);
        I2 = I2_p(s+1:end-s,s+1:end-s,:);
    else
        I1 = I1_p(s+1:end-s,s+1:end-s);
        I2 = I2_p(s+1:end-s,s+1:end-s);
    end
    loc1 = loc1-s;
    loc2 = loc2-s;
    
    im3 = rr_appendimages(I1,I2);
    % h = figure('visible',visible);
    if(nargin>6)
        h = figure('NumberTitle', 'off', 'Name', titles,'visible',visible);
    else
        h = figure('visible',visible);
    end
    imshow(im3,[])
    title(['matched points (' num2str(size(loc1,1)) ')']);
    hold on
    cols = size(I1,2);
    
    for i=1:size(loc1,1)
        line([loc1(i,2) loc2(i,2)+cols],[loc1(i,1) loc2(i,1)], 'Color', 'y');
        plot(loc1(i,2),loc1(i,1),'g.')
        plot(loc1(i,2),loc1(i,1),'go')
        plot(loc2(i,2)+cols,loc2(i,1),'g.')
        plot(loc2(i,2)+cols,loc2(i,1),'go')
    end
    hold off
    
    %     saveas(h,['E:\学习\ur-sift实验\ur-sift\paper\fig6 组图\' num2str((imageNum*2-1)) '-' num2str((imageNum*2)) '.bmp']);
end
end