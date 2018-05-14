%多元均值向量
Mean (magic04)
%the inner-product(内积计算）
 Z=magic04-ones*means

X=Z'*Z*(1/19025)
%the outer-product（外部乘积）
sum=0
for i=1:19025
S=Z(i,:)'*Z(i,:)
sum=sum+S*(1/19025)
end

%求夹角余弦 
 A=magic04(:,1) %提取特征向量
B=magic04(:,2)
c=dot(A,B)/(norm(A)*norm(B))

%散点图 
x=1:length(A);
plot(x,A,'go','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
hold on
plot(x,B,'cs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','c','MarkerSize',10)
hold on
axis([0 100 0 100])
grid on
T=legend('1#','2#');
set(T,'Fontsize',12);

%假设属性1是正态分布的，绘制其概率密度函数。 

%均值：
c=mean(A)
%标准差：
sigma=std(A)
  x = -200:0.1:200;
y = normpdf(x, c, sigma);
plot(x,y);
grid on;

%鸢尾花 
%线性核函数 the linear kernel ：
 for i=1:4
for j=1:4
k(i,j)=iris(:,i)'*iris(:,j)
end
end
%齐次二次核the homogeneous quadratic kernel ：
 for i=1:4
for j=1:4
k(i,j)=(iris(:,i)'*iris(:,j))^2
end
 end

%the linear kernel进行中心化：
 Y=eye(4,4)-eye(4,4)/4%生成对角线为1的4*4 矩阵
 centered=Y*k*Y %中心化
 
 
 %对the homogeneous quadratic kernel中心化
 centered2=Y*(k^2)*Y

 
% 对the linear kernel进行归一化
W=diag(diag(k))%取对角线元素并生成矩阵
normalized=W^(1/2)*k*W^(1/2)%归一化

%对the homogeneous quadratic kernel归一化：
normalized=W^(1/2)*(k^2)*W^(1/2)

%每个点x变换到特征空间?（x）
%取 iris的每一列数据 
X1=iris(:,1)
X2=iris(:,2)
X3=iris(:,3)
x4=iris(:,4)
%变化到特征空间?（x）
X1=x1/sqrt(sum(x1.*x1))
X2=x2/sqrt(sum(x2.*x2))
X3=x3/sqrt(sum(x3.*x3))
X4=x4/sqrt(sum(x4.*x4))
%合并矩阵
 iris2=[x1,x2,x3,x4]

 for i=1:4%线性核函数
for j=1:4
k2(i,j)=iris2(:,i)'*iris2(:,j)
end
 end

 %中心化 ：
 centered3=Y*k2*Y
 
 %归一化
 
  W=diag(diag(k2))
   normalized2=W^(1/2)*k2*W^(1/2)


