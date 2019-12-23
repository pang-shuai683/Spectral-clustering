%%%%%%%%%%%%%% �׾��� K-means
% load 'data.mat'
% d1=allpts(1:20,:)';         %��ǩΪ1������
% d2=allpts(21:50,:)';        %��ǩΪ2������
% %%%%%%%%%%%%%%%%%%%%
% 
% cluster1=d1';
% cluster2=d2';
data = xlsread( 'C:\Users\cloud\Desktop\��ҵ����\ʵ�����\SSR����汾\data\AEEEM\PDE.xlsx');
[N,m] = size(data);
data1=data(:,2:m-1);

%����Ԥ���� ��һ��
Xsq = data1.*data1;      
divmat=repmat(sqrt(sum(Xsq')'),1,m-2);
allpts=data1./divmat;

%% compute A (step 1)
%%%%  experiment with sigsq in question 8
sigsq=.9;    %�������ƶȾ���W������Ǹ�����
Aisq=sum(allpts.*allpts,2);
Dotprod=allpts*allpts';
distmat=-repmat(Aisq',N,1)-repmat(Aisq,1,N)+2*Dotprod;
Afast=exp(distmat/(2*sigsq));   %���ƶȾ���
A = Afast-diag(diag(Afast));     %���Խ�����Ϊ��


%% step 2
D = diag(sum(A')) ;  %�Ⱦ���
L=D^(-.5)*A*D^(-.5) ;  %������˹����

%%������ֵ����������������ǰK��
K = 10;   %��ά ������KΪ2
[X,di]=eig(L);
[Xsort,Dsort]=eigsort(X,di);
Xuse=Xsort(:,1:K);   

%%�Եõ���K�����������ľ�����й�һ�� ��һ�����ǽ������������ȱ��1
Xsq = Xuse.*Xuse;      
divmat=repmat(sqrt(sum(Xsq')'),1,K);
Y=Xuse./divmat;

%% ����K-means����
[c,Dsum,z] = kmeans(Y,2);  

kk=c;
c1=find(kk==1);
c2=find(kk==2);