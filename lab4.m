X=iris(:,3:4);  
  
%%KNN k distance graph, to determine the epsilon  
A=X;  
numData=size(A,1);  
Kdist=zeros(numData,1);  
[IDX,Dist]=knnsearch(A(2:numData,:),A(1,:));  
Kdist(1)=Dist;  
for i=2:size(A,1)  
    [IDX,Dist] = knnsearch(A([1:i-1,i+1:numData],:),A(i,:));  
    Kdist(i)=Dist;  
end  
[sortKdist,sortKdistIdx]=sort(Kdist,'descend');  
distX=[1:numData]';  
plot(distX,sortKdist,'r+-','LineWidth',2);  
set(gcf,'position',[1000 340 350 350]);  
grid on;  
  
%% Run DBSCAN Clustering Algorithm  
epsilon= 0.15 ;  
MinPts=  3   ;  
IDX1=DBSCAN(X,epsilon,MinPts);  
%% Plot Results  
figure;  
PlotClusterinResult(X, IDX1);  
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);  
set(gcf,'position',[30 -10 500 500]);   
  
  
epsilon= 0.25 ;  
MinPts=  3   ;  
IDX2=DBSCAN(X,epsilon,MinPts);  
%% Plot Results  
figure;  
PlotClusterinResult(X, IDX2);  
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);  
set(gcf,'position',[530 -10 500 500]);  
  
epsilon= 0.5 ;  
MinPts=  3   ;  
IDX3=DBSCAN(X,epsilon,MinPts);  
%% Plot Results  
figure;  
PlotClusterinResult(X, IDX3);  
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);  
set(gcf,'position',[30 380 500 500]);  
  
  
%DBSCAN算法子函数，需另外创建.m文件保存  
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)  
% All rights reserved. Please read the "license.txt" for license terms.  
%  
% Project Code: YPML110  
% Project Title: Implementation of DBSCAN Clustering in MATLAB  
% Publisher: Yarpiz (www.yarpiz.com)  
%   
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)  
%   
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com  
function [IDX, isnoise]=DBSCAN(X,epsilon,MinPts)  
    C=0;  
    n=size(X,1);  
    IDX=zeros(n,1);  
    D=pdist2(X,X);  
    visited=false(n,1);  
    isnoise=false(n,1);  
    for i=1:n  
        if ~visited(i)  
            visited(i)=true;  
              
            Neighbors=RegionQuery(i);  
            if numel(Neighbors)<MinPts  
                % X(i,:) is NOISE  
                isnoise(i)=true;  
            else  
                C=C+1;  
                ExpandCluster(i,Neighbors,C);  
            end  
              
        end  
    end  
      
    function ExpandCluster(i,Neighbors,C)  
        IDX(i)=C;  
          
        k = 1;  
        while true  
            j = Neighbors(k);  
              
            if ~visited(j)  
                visited(j)=true;  
                Neighbors2=RegionQuery(j);  
                if numel(Neighbors2)>=MinPts  
                    Neighbors=[Neighbors Neighbors2];   %#ok  
                end  
            end  
            if IDX(j)==0  
                IDX(j)=C;  
            end  
              
            k = k + 1;  
            if k > numel(Neighbors)  
                break;  
            end  
        end  
    end  
      
    function Neighbors=RegionQuery(i)  
        Neighbors=find(D(i,:)<=epsilon);  
    end  
  
end  
  
%结果显示子函数，需另外创建.m文件保存  
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)  
% All rights reserved. Please read the "license.txt" for license terms.  
%  
% Project Code: YPML110  
% Project Title: Implementation of DBSCAN Clustering in MATLAB  
% Publisher: Yarpiz (www.yarpiz.com)  
%   
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)  
%   
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com  
  
function PlotClusterinResult(X, IDX)  
  
    k=max(IDX);  
  
    Colors=hsv(k);  
  
    Legends = {};  
    for i=0:k  
        Xi=X(IDX==i,:);  
        if i~=0  
            Style = 'x';  
            MarkerSize = 8;  
            Color = Colors(i,:);  
            Legends{end+1} = ['Cluster #' num2str(i)];  
        else  
            Style = 'o';  
            MarkerSize = 6;  
            Color = [0 0 0];  
            if ~isempty(Xi)  
                Legends{end+1} = 'Noise';  
            end  
        end  
        if ~isempty(Xi)  
            plot(Xi(:,1),Xi(:,2),Style,'MarkerSize',MarkerSize,'Color',Color);  
        end  
        hold on;  
    end  
    hold off;  
    axis equal;  
    grid on;  
    legend(Legends);  
    legend('Location', 'NorthEastOutside');  
end 