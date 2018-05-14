 data=iris;
% ���ɾ�����  
createTree(data);
///////////
%creatTree
function [output_args] = createTree( data )  
    [m,n] = size(data);  
    disp('original data:');  
    disp(data);  
    classList = data(:,n);  
    classOne = 1;%��¼��һ����ĸ���  
    for i = 2:m  
        if classList(i,:) == classList(1,:)  
            classOne = classOne+1;  
        end  
    end  
      
    % ���ȫ��ͬ  
    if classOne == m  
        disp('final data: ');  
        disp(data);  
        return;  
    end  
      
    % ����ȫ������  
    if n == 1  
        disp('final data: ');  
        disp(data);  
        return;  
    end  
      
    bestFeat = chooseBestFeature(data);  
    disp(['bestFeat: ', num2str(bestFeat)]);  
    featValues = unique(data(:,bestFeat));  
    numOfFeatValue = length(featValues);  
      
    for i = 1:numOfFeatValue  
        createTree(splitData(data, bestFeat, featValues(i,:)));  
        disp('-------------------------');  
    end  
end  
//////////
%% ѡ����Ϣ������������  
function [ bestFeature ] = chooseBestFeature( data )  
    [m,n] = size(data);% �õ����ݼ��Ĵ�С  
      
    % ͳ�������ĸ���  
    numOfFeatures = n-1;%���һ�������  
    % ԭʼ����  
    baseEntropy = calEntropy(data);  
      
    bestInfoGain = 0;%��ʼ����Ϣ����  
    bestFeature = 0;% ��ʼ����ѵ�����λ  
      
    % ��ѡ��ѵ�����λ  
    for j = 1:numOfFeatures  
        featureTemp = unique(data(:,j));  
        numF = length(featureTemp);%���Եĸ���  
        newEntropy = 0;%����֮�����  
        for i = 1:numF  
            subSet = splitData(data, j, featureTemp(i,:));  
            [m_1, n_1] = size(subSet);  
            prob = m_1./m;  
            newEntropy = newEntropy + prob * calEntropy(subSet);  
        end  
          
        %��������  
        infoGain = baseEntropy - newEntropy;  
          
        if infoGain > bestInfoGain  
            bestInfoGain = infoGain;  
            bestFeature = j;  
        end  
    end  
end  
//////////
function [ entropy ] = calEntropy( data )  
    [m,n] = size(data);  
      
    % �õ�������  
    label = data(:,n);  
      
    % �������label  
    label_deal = unique(label);  
      
    numLabel = length(label_deal);  
    prob = zeros(numLabel,2);  
      
    % ͳ�Ʊ�ǩ  
    for i = 1:numLabel  
        prob(i,1) = label_deal(i,:);  
        for j = 1:m  
            if label(j,:) == label_deal(i,:)  
                prob(i,2) = prob(i,2)+1;  
            end  
        end  
    end  
      
    % ������  
    prob(:,2) = prob(:,2)./m;  
    entropy = 0;  
    for i = 1:numLabel  
        entropy = entropy - prob(i,2) * log2(prob(i,2));  
    end  
end  

//////
function [ subSet ] = splitData( data, axis, value )  
    [m,n] = size(data);%�õ����������ݵĴ�С  
      
    subSet = data;  
    subSet(:,axis) = [];  
    k = 0;  
    for i = 1:m  
        if data(i,axis) ~= value  
            subSet(i-k,:) = [];  
            k = k+1;  
        end  
    end     
end  