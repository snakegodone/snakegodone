%��Ԫ��ֵ����
Mean (magic04)
%the inner-product(�ڻ����㣩
 Z=magic04-ones*means

X=Z'*Z*(1/19025)
%the outer-product���ⲿ�˻���
sum=0
for i=1:19025
S=Z(i,:)'*Z(i,:)
sum=sum+S*(1/19025)
end

%��н����� 
 A=magic04(:,1) %��ȡ��������
B=magic04(:,2)
c=dot(A,B)/(norm(A)*norm(B))

%ɢ��ͼ 
x=1:length(A);
plot(x,A,'go','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
hold on
plot(x,B,'cs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','c','MarkerSize',10)
hold on
axis([0 100 0 100])
grid on
T=legend('1#','2#');
set(T,'Fontsize',12);

%��������1����̬�ֲ��ģ�����������ܶȺ����� 

%��ֵ��
c=mean(A)
%��׼�
sigma=std(A)
  x = -200:0.1:200;
y = normpdf(x, c, sigma);
plot(x,y);
grid on;

%�β�� 
%���Ժ˺��� the linear kernel ��
 for i=1:4
for j=1:4
k(i,j)=iris(:,i)'*iris(:,j)
end
end
%��ζ��κ�the homogeneous quadratic kernel ��
 for i=1:4
for j=1:4
k(i,j)=(iris(:,i)'*iris(:,j))^2
end
 end

%the linear kernel�������Ļ���
 Y=eye(4,4)-eye(4,4)/4%���ɶԽ���Ϊ1��4*4 ����
 centered=Y*k*Y %���Ļ�
 
 
 %��the homogeneous quadratic kernel���Ļ�
 centered2=Y*(k^2)*Y

 
% ��the linear kernel���й�һ��
W=diag(diag(k))%ȡ�Խ���Ԫ�ز����ɾ���
normalized=W^(1/2)*k*W^(1/2)%��һ��

%��the homogeneous quadratic kernel��һ����
normalized=W^(1/2)*(k^2)*W^(1/2)

%ÿ����x�任�������ռ�?��x��
%ȡ iris��ÿһ������ 
X1=iris(:,1)
X2=iris(:,2)
X3=iris(:,3)
x4=iris(:,4)
%�仯�������ռ�?��x��
X1=x1/sqrt(sum(x1.*x1))
X2=x2/sqrt(sum(x2.*x2))
X3=x3/sqrt(sum(x3.*x3))
X4=x4/sqrt(sum(x4.*x4))
%�ϲ�����
 iris2=[x1,x2,x3,x4]

 for i=1:4%���Ժ˺���
for j=1:4
k2(i,j)=iris2(:,i)'*iris2(:,j)
end
 end

 %���Ļ� ��
 centered3=Y*k2*Y
 
 %��һ��
 
  W=diag(diag(k2))
   normalized2=W^(1/2)*k2*W^(1/2)


