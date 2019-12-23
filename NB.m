function [performance_all] = NB(predictlabel, targetlabel,loc)


targetlabel = targetlabel-1;
% sourcedata =source(:,1:end-1);
% sourcelabel = source(:,end);
% targetdata = target(:,1:end-1);
% targetlabel = target(:,end);
% 
predictlabel = predictlabel';
predictlabel = predictlabel-1;
% [predictlabel targetlabel]
% 
% Factor = NaiveBayes.fit(sourcedata, sourcelabel,'Distribution','normal');
% [Score,predictlabel] = posterior(Factor, targetdata);
% 
% score =Score(:,2);
% logitFit = score;

predictlabel_zhuanzhi = predictlabel';
targetlabel_zhuanzhi = targetlabel';

TN = 0;
FP = 0;
FN = 0;
TP = 0;
for i = 1:size(targetlabel_zhuanzhi,2)
    if (targetlabel_zhuanzhi(i)== 0 && predictlabel_zhuanzhi(i) == 0)
	    TN = TN + 1;
	elseif (targetlabel_zhuanzhi(i) == 0 && predictlabel_zhuanzhi(i) == 1)
        FP = FP + 1;
	elseif (targetlabel_zhuanzhi(i) == 1 && predictlabel_zhuanzhi(i) == 0)
        FN = FN + 1;
    else
        TP = TP + 1;
    end
end
r = [TN , FP , FN , TP];
%�����ָ��
Accuracy = (TN + TP) / (TN + FP + FN + TP);
Recall = TP / (TP + FN);
Pd = Recall;
Precision = TP / (TP + FP);
False_Positive = FP / (FP + TN);
Pf = False_Positive;
F_measure = 2 * Precision * Recall / (Precision + Recall);
F_2 = 5 * Precision * Recall / (4 * Precision + Recall);
G_measure = 2 * Recall * (1 - False_Positive) / (Recall + (1 - False_Positive));
g_mean = sqrt((TN / (TN + FP)) * (TP / (TP + FN)));
Bal = 1- sqrt((0-Pf)^2+(1-Pd)^2)/sqrt(2);
MCC = (TP * TN - FN * FP) / sqrt((TP + FN) * (TP + FP) * (FN + TN) * (FP + TN));
%  [~,~,~,AUC] = perfcurve(targetlabel,logitFit,1);

% Performance_tradition = [Pd, Pf, Precision, Recall, F_measure, F_2, G_measure, g_mean, Bal, MCC, AUC];
Performance_tradition = [F_measure, F_2, G_measure, g_mean, Bal, MCC];

%����effort-aware������ָ��
index_defect = find(predictlabel == 1); % Ԥ��Ϊ�д��ģ�����
index_non_defect = find(predictlabel == 0);
predicted_defect_actuallabel = targetlabel(index_defect); % Ԥ��Ϊ��ȱ�ݵ�ģ�����ʵlabel
predicted_non_defect_actuallabel = targetlabel(index_non_defect);% ��Ԥ��Ϊû��ȱ�ݵ�ģ�����ʵlabel
predicted_defect_LOC = loc(index_defect); % Ԥ��Ϊ��ȱ��ģ���LOCֵ
predicted_non_defect_LOC = loc(index_non_defect);% Ԥ��Ϊû��ȱ��ģ���LOCֵ

[LOC_defect_SORT, index_defect_LOC] = sort(predicted_defect_LOC); % Ԥ��Ϊ��ȱ�ݵ�ģ�鰴��LOC������������
[LOC_non_defect_SORT, index_non_defect_LOC] = sort(predicted_non_defect_LOC); % Ԥ��Ϊû��ȱ�ݵ�ģ�鰴��LOC������������
LOC_SORT = [LOC_defect_SORT;LOC_non_defect_SORT]; % Ԥ��Ϊ��ȱ�ݵĺ�Ԥ��Ϊû��ȱ�ݵİ���LOC����ƴ������

predicted_defect_actuallabel_sort = predicted_defect_actuallabel(index_defect_LOC); %Ԥ��Ϊ��ȱ�ݵ�ģ�����ʵ��ǩ����LOC����
predicted_non_defect_actuallabel_sort = predicted_non_defect_actuallabel(index_non_defect_LOC);%Ԥ��Ϊû��ȱ�ݵ�ģ�����ʵ��ǩ����LOC����
predicted_actuallabel_sort = [predicted_defect_actuallabel_sort;predicted_non_defect_actuallabel_sort]; %ƴ��

CumLOC = 100 * cumsum(LOC_SORT)./sum(loc); %�����ۻ����Ĵ������
Cutoffindex = find(CumLOC >= 20);          %�����20%�Ĵ���ʱֹͣ

% Cutoffindex
Cutoffindex(1);

predicted_actuallabel_cutoff20 = predicted_actuallabel_sort(1:Cutoffindex(1));

% Suppose having a dataset with M modules and N defects.
% After inspecting 20% LOC, we inspected m modules and found n defects. 
% Besides, when we find the first defective module, we have inspected k modules.

M = size(targetlabel,1);
N = sum(targetlabel);
m = Cutoffindex(1);
n = sum(predicted_actuallabel_cutoff20);
k_index = find(predicted_actuallabel_cutoff20 == 1);
if n ~= 0
    k = k_index(1);
else
    k = NaN;
end

% ����ָ��
Effort_Recall = n/N;
Effort_Precision = n/m;
Effort_F = 2 * Effort_Precision * Effort_Recall / (Effort_Precision + Effort_Recall);
Effort_F_2 = 5 * Effort_Precision * Effort_Recall / (4 * Effort_Precision + Effort_Recall);
PCI20 = m/M;
IFA = k;

% Performance_novel = [Effort_Recall, Effort_Precision, Effort_F, Effort_F_2, PCI20, IFA];
Performance_EA = [Effort_Precision, Effort_Recall, Effort_F, Effort_F_2, PCI20, IFA];

performance_all = [Performance_tradition, Performance_EA];