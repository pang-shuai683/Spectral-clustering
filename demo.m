%%%%%%%%%%%%%% 谱聚类 K-means
% load 'data.mat'
% d1=allpts(1:20,:)';         %标签为1的数据
% d2=allpts(21:50,:)';        %标签为2的数据
% %%%%%%%%%%%%%%%%%%%%
% 
% cluster1=d1';
% cluster2=d2';
data = xlsread( 'C:\Users\cloud\Desktop\毕业论文\实验代码\SSR会议版本\data\AEEEM\PDE.xlsx');
[N,m] = size(data);
data1=data(:,2:m-1);

%数据预处理 归一化
Xsq = data1.*data1;      
divmat=repmat(sqrt(sum(Xsq')'),1,m-2);
allpts=data1./divmat;

%% compute A (step 1)
%%%%  experiment with sigsq in question 8
sigsq=.9;    %计算相似度矩阵W里面的那个参数
Aisq=sum(allpts.*allpts,2);
Dotprod=allpts*allpts';
distmat=-repmat(Aisq',N,1)-repmat(Aisq,1,N)+2*Dotprod;
Afast=exp(distmat/(2*sigsq));   %相似度矩阵
A = Afast-diag(diag(Afast));     %将对角线设为零


%% step 2
D = diag(sum(A')) ;  %度矩阵
L=D^(-.5)*A*D^(-.5) ;  %拉普拉斯矩阵

%%找特征值特征向量并排序，找前K个
K = 10;   %降维 特征数K为2
[X,di]=eig(L);
[Xsort,Dsort]=eigsort(X,di);
Xuse=Xsort(:,1:K);   

%%对得到的K个特征向量的矩阵进行归一化 归一化就是将特征向量长度变成1
Xsq = Xuse.*Xuse;      
divmat=repmat(sqrt(sum(Xsq')'),1,K);
Y=Xuse./divmat;

%% 进行K-means聚类
[c,Dsum,z] = kmeans(Y,2);  

kk=c;
c1=find(kk==1);
c2=find(kk==2);