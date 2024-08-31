clear
clc

load('MtxGroup1_mother_with_demographic_data.mat','MtxGroup1_mother');


%%

 M=MtxGroup1_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));    
    ix=find((M.SP_Ct_Mean(indx,:))>0);    
    ix2=find((M.RSV_CT(indx,:))<99);    
    mean_SP(j)= mean( (M.SP_Ct_Mean((indx(ix)),:) ) );
    sd_SP(j)=  std( (M.SP_Ct_Mean((indx(ix)),:) ) );    
    mean_RSV(j)= mean( (M.RSV_CT((indx(ix2)),:) ) );
    %     sd_RSV(j)=  std( (M.RSV_CT((indx(ix2)),:) ) );
    mother_HIV(j)=M.Mother_HIV_status((indx(ix2)),:) ;
    birth_weight(j)=M.Birth_Weight((indx(ix2)),:) ;
    mother_age(j)=M.Mother_age((indx(ix2)),:) ;
    age_in_days_acquisition_SP(j)= time2num((M.Date_of_Visit(indx(ix(1))) - M.Child_DOB(indx(ix(1)))),"days");    
end

group=ones(1,length(subjects))*1;

Mt1=table(subjects',mean_SP',sd_SP',mean_RSV',mother_HIV',birth_weight',mother_age',group',age_in_days_acquisition_SP');
Mt1=renamevars(Mt1,"Var1","Subjects");
Mt1=renamevars(Mt1,"Var2","mean_SP");
Mt1=renamevars(Mt1,"Var3","sd_SP");
Mt1=renamevars(Mt1,"Var4","mean_RSV");
Mt1=renamevars(Mt1,"Var5","mother_HIV");
Mt1=renamevars(Mt1,"Var6","birth_weight");
Mt1=renamevars(Mt1,"Var7","mother_age");
Mt1=renamevars(Mt1,"Var8","group");
Mt1=renamevars(Mt1,"Var9","age_days_SP_acquisition");

%%
clearvars -except Mt1

load('MtxGroup2_mother_with_demographic_data.mat','MtxGroup2_mother');


M=MtxGroup2_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));    
    ix=find((M.SP_Ct_Mean(indx,:))>0);    
    ix2=find((M.RSV_CT(indx,:))<99);    
    mean_SP(j)= mean( (M.SP_Ct_Mean((indx(ix)),:) ) );
    sd_SP(j)=  std( (M.SP_Ct_Mean((indx(ix)),:) ) );    
    mean_RSV(j)= mean( (M.RSV_CT((indx(ix2)),:) ) );
    %     sd_RSV(j)=  std( (M.RSV_CT((indx(ix2)),:) ) );
    mother_HIV(j)=M.Mother_HIV_status((indx(ix2)),:) ;
    birth_weight(j)=M.Birth_Weight((indx(ix2)),:) ;
    mother_age(j)=M.Mother_age((indx(ix2)),:) ;    
    age_in_days_acquisition_SP(j)= time2num((M.Date_of_Visit(indx(ix(1))) - M.Child_DOB(indx(ix(1)))),"days");
end

group=ones(1,length(subjects))*2;

Mt2=table(subjects',mean_SP',sd_SP',mean_RSV',mother_HIV',birth_weight',mother_age',group',age_in_days_acquisition_SP');
Mt2=renamevars(Mt2,"Var1","Subjects");
Mt2=renamevars(Mt2,"Var2","mean_SP");
Mt2=renamevars(Mt2,"Var3","sd_SP");
Mt2=renamevars(Mt2,"Var4","mean_RSV");
Mt2=renamevars(Mt2,"Var5","mother_HIV");
Mt2=renamevars(Mt2,"Var6","birth_weight");
Mt2=renamevars(Mt2,"Var7","mother_age");
Mt2=renamevars(Mt2,"Var8","group");
Mt2=renamevars(Mt2,"Var9","age_days_SP_acquisition");

%%

clearvars -except Mt1 Mt2

load('MtxGroup3_mother_with_demographic_data.mat','MtxGroup3_mother');


M=MtxGroup3_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));
    
     ix=find((M.SP_Ct_Mean(indx,:))>0);
     
     ix2=find((M.RSV_CT(indx,:))<99);
    
    mean_SP(j)= mean( (M.SP_Ct_Mean((indx(ix)),:) ) );
    sd_SP(j)=  std( (M.SP_Ct_Mean((indx(ix)),:) ) );
    
    mean_RSV(j)= mean( (M.RSV_CT((indx(ix2)),:) ) );
%     sd_RSV(j)=  std( (M.RSV_CT((indx(ix2)),:) ) );
    mother_HIV(j)=M.Mother_HIV_status((indx(1)),:) ;
    birth_weight(j)=M.Birth_Weight((indx(1)),:) ;
    mother_age(j)=M.Mother_age((indx(1)),:) ;
    
    age_in_days_acquisition_SP(j)= time2num((M.Date_of_Visit(indx(ix(1))) - M.Child_DOB(indx(ix(1)))),"days");
    
    
    
end

group=ones(1,length(subjects))*3;

Mt3=table(subjects',mean_SP',sd_SP',mean_RSV',mother_HIV',birth_weight',mother_age',group',age_in_days_acquisition_SP');
Mt3=renamevars(Mt3,"Var1","Subjects");
Mt3=renamevars(Mt3,"Var2","mean_SP");
Mt3=renamevars(Mt3,"Var3","sd_SP");
Mt3=renamevars(Mt3,"Var4","mean_RSV");
Mt3=renamevars(Mt3,"Var5","mother_HIV");
Mt3=renamevars(Mt3,"Var6","birth_weight");
Mt3=renamevars(Mt3,"Var7","mother_age");
Mt3=renamevars(Mt3,"Var8","group");
Mt3=renamevars(Mt3,"Var9","age_days_SP_acquisition");

%%


clearvars -except Mt1 Mt2 Mt3 

load('MtxGroup4_mother_with_demographic_data.mat','MtxGroup4_mother');


M=MtxGroup4_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));
    
     ix=find((M.SP_Ct_Mean(indx,:))>0);
     
     ix2=find((M.RSV_CT(indx,:))<99);
    
    mean_SP(j)= mean( (M.SP_Ct_Mean((indx(ix)),:) ) );
    sd_SP(j)=  std( (M.SP_Ct_Mean((indx(ix)),:) ) );
    
    mean_RSV(j)= mean( (M.RSV_CT((indx(ix2)),:) ) );

    mother_HIV(j)=M.Mother_HIV_status(j,:) ;
    birth_weight(j)=M.Birth_Weight(j,:) ;
    
    mother_age(j)=M.Mother_age((indx(1)),:) ;
end

group=ones(1,length(subjects))*5;

% mean_RSV=zeros(1,length(subjects));


Mt4=table(subjects',mean_SP',sd_SP',mean_RSV',mother_HIV',birth_weight',mother_age',group');
Mt4=renamevars(Mt4,"Var1","Subjects");
Mt4=renamevars(Mt4,"Var2","mean_SP");
Mt4=renamevars(Mt4,"Var3","sd_SP");
Mt4=renamevars(Mt4,"Var4","mean_RSV");
Mt4=renamevars(Mt4,"Var5","mother_HIV");
Mt4=renamevars(Mt4,"Var6","birth_weight");
Mt4=renamevars(Mt4,"Var7","mother_age");
Mt4=renamevars(Mt4,"Var8","group");




%%

clearvars -except Mt1 Mt2 Mt3 Mt4

load('MtxGroup5_mother_with_demographic_data.mat','MtxGroup5_mother');


M=MtxGroup5_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));
    
     ix=find((M.SP_Ct_Mean(indx,:))>0);
     
     ix2=find((M.RSV_CT(indx,:))==99);
    
    mean_SP(j)= mean( (M.SP_Ct_Mean((indx(ix)),:) ) );
    sd_SP(j)=  std( (M.SP_Ct_Mean((indx(ix)),:) ) );
    
%     mean_RSV(j)= mean( (M.RSV_CT((indx(ix)),:) ) );

    mother_HIV(j)=M.Mother_HIV_status(j,:) ;
    birth_weight(j)=M.Birth_Weight(j,:) ;
    
    mother_age(j)=M.Mother_age((indx(1)),:) ;
    
    age_in_days_acquisition_SP(j)= time2num((M.Date_of_Visit(indx(ix(1))) - M.Child_DOB(indx(ix(1)))),"days");

end

group=ones(1,length(subjects))*5;

mean_RSV=zeros(1,length(subjects));


Mt5=table(subjects',mean_SP',sd_SP',mean_RSV',mother_HIV',birth_weight',mother_age',group',age_in_days_acquisition_SP');
Mt5=renamevars(Mt5,"Var1","Subjects");
Mt5=renamevars(Mt5,"Var2","mean_SP");
Mt5=renamevars(Mt5,"Var3","sd_SP");
Mt5=renamevars(Mt5,"Var4","mean_RSV");
Mt5=renamevars(Mt5,"Var5","mother_HIV");
Mt5=renamevars(Mt5,"Var6","birth_weight");
Mt5=renamevars(Mt5,"Var7","mother_age");
Mt5=renamevars(Mt5,"Var8","group");
Mt5=renamevars(Mt5,"Var9","age_days_SP_acquisition");



%%

clearvars -except Mt1 Mt2 Mt3 Mt4 Mt5

load('MtxGroup6_mother_with_demographic_data.mat','MtxGroup6_mother');


M=MtxGroup6_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));
    
     ix=find((M.SP_Ct_Mean(indx,:))>0);
     
     ix2=find((M.RSV_CT(indx,:))==99);
    
    mean_SP(j)= mean( (M.SP_Ct_Mean((indx(ix)),:) ) );
    sd_SP(j)=  std( (M.SP_Ct_Mean((indx(ix)),:) ) );
    
%     mean_RSV(j)= mean( (M.RSV_CT((indx(ix)),:) ) );

    mother_HIV(j)=M.Mother_HIV_status(j,:) ;
    birth_weight(j)=M.Birth_Weight(j,:) ;
    
    mother_age(j)=M.Mother_age((indx(1)),:) ;
end

group=ones(1,length(subjects))*6;

mean_RSV=zeros(1,length(subjects));


Mt6=table(subjects',mean_SP',sd_SP',mean_RSV',mother_HIV',birth_weight',mother_age',group');
Mt6=renamevars(Mt6,"Var1","Subjects");
Mt6=renamevars(Mt6,"Var2","mean_SP");
Mt6=renamevars(Mt6,"Var3","sd_SP");
Mt6=renamevars(Mt6,"Var4","mean_RSV");
Mt6=renamevars(Mt6,"Var5","mother_HIV");
Mt6=renamevars(Mt6,"Var6","birth_weight");
Mt6=renamevars(Mt6,"Var7","mother_age");
Mt6=renamevars(Mt6,"Var8","group");


%%
col_bl=[0 0.4470 0.7410];
col_org=[0.8500 0.3250 0.0980];
col_yl=[0.9290 0.6940 0.1250];
col_prpl=[0.4940 0.1840 0.5560];
col_blk =[0 0 0];

% n=max([length(Mt1.age_days_SP_acquisition),length(Mt2.age_days_SP_acquisition),length(Mt3.age_days_SP_acquisition),length(Mt5.age_days_SP_acquisition)]);
% Mt1.age_days_SP_acquisition(end+1:n)=nan;
% Mt2.age_days_SP_acquisition(end+1:n)=nan;
% Mt3.age_days_SP_acquisition(end+1:n)=nan;
% Mt5.age_days_SP_acquisition(end+1:n)=nan;


X=[Mt1.age_days_SP_acquisition;Mt2.age_days_SP_acquisition;Mt3.age_days_SP_acquisition;Mt5.age_days_SP_acquisition];
XX=time2num(X,"days");
Y=[ones(length(Mt1.age_days_SP_acquisition),1);2*ones(length(Mt2.age_days_SP_acquisition),1);3*ones(length(Mt3.age_days_SP_acquisition),1);4*ones(length(Mt5.age_days_SP_acquisition),1)];
figure
xSize = 15; Xs=xSize; ySize = 11.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

% histogram(Mt5.age_days_SP_acquisition,30)
boxplot(XX,Y,'Color','k','Symbol','');

hold on
scatter(ones(length(Mt1.age_days_SP_acquisition),1),Mt1.age_days_SP_acquisition,1000,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(length(Mt2.age_days_SP_acquisition),1),Mt2.age_days_SP_acquisition,1000,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(3*ones(length(Mt3.age_days_SP_acquisition),1),Mt3.age_days_SP_acquisition,1000,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(4*ones(length(Mt5.age_days_SP_acquisition),1),Mt5.age_days_SP_acquisition,1000,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot(1,mean(Mt1.age_days_SP_acquisition),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
plot(2,mean(Mt2.age_days_SP_acquisition),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
plot(3,mean(Mt3.age_days_SP_acquisition),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
plot(4,mean(Mt5.age_days_SP_acquisition),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
hold on
xlim([0.5,4.5])
set(gca,'Fontsize',40);box on;
xticks((1:1:4))
set(gca, 'XTickLabel', {'SP$\rightarrow$RSV','SP\&RSV','RSV$\rightarrow$SP','SP'},'TickLabelInterpreter','latex','Fontsize',25)
ylim([-5,119])
ylabel('SP acquisition (days)','Interpreter','latex')
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(30)

% % plots line to denote significant difference
% yt = get(gca, 'YTick');
% axis([xlim    0  ceil(max(yt)*1.1)])
% % xt = get(gca, 'XTick');
% xt=[1,2,3,4];
% 
% hold on
% plot(xt([2 3]), [1 1]*max(yt)*1.10, 'k-','LineWidth',6)
% hold on
% plot(xt([3 4]), [1 1]*max(yt)*1.19, 'k-','LineWidth',6)
% hold on
% plot(xt([1 3]), [1 1]*max(yt)*1.15, 'k-','LineWidth',6)
% 
% dim1 = [.51 .6 .9 .3];
% dim2 = [.7 .65 .9 .3];
% dim3 = [.42 .63 .9 .3];
% 
% 
% name={'$ * $'};
% annotation('textbox',dim1,'String',name,'interpreter','latex','Fontsize',55,'Color', 'k','EdgeColor','none');
% hold on
% annotation('textbox',dim2,'String',name,'interpreter','latex','Fontsize',55,'Color', 'k','EdgeColor','none');
% hold on
% annotation('textbox',dim3,'String',name,'interpreter','latex','Fontsize',55,'Color', 'k','EdgeColor','none');

ylim([-5,125])
set(gca,'Fontsize',60);box on;

saveas(gcf,'boxplot_SP_acquisition_mothers_v2.png');

%% test for normality

[hc1,pc1] = kstest(Mt1.age_days_SP_acquisition);
[hc2,pc2] = kstest(Mt2.age_days_SP_acquisition);
[hc3,pc3] = kstest(Mt3.age_days_SP_acquisition);
[hc4,pc4] = kstest(Mt5.age_days_SP_acquisition);

pc=[pc1,pc2,pc3,pc4];

% Check for significant differences
alpha = 0.05;  % Set significance level
for i=1:4
    if pc(i) < alpha
        fprintf('Infants Group %d is NOT normally distributed.\n', i );
    else
        fprintf('Infants Group %d is normally distributed.\n', i );
    end
end


%% T-test NOT recommended as samples are NOT normally distributed

% [h1,p1]=ttest2(Mt1.age_days_SP_acquisition,Mt2.age_days_SP_acquisition);
% [h2,p2]=ttest2(Mt1.age_days_SP_acquisition,Mt3.age_days_SP_acquisition); %%%%%%
% [h3,p3]=ttest2(Mt1.age_days_SP_acquisition,Mt5.age_days_SP_acquisition);
% [h1a,p1a]=ttest2(Mt2.age_days_SP_acquisition,Mt3.age_days_SP_acquisition); %%%%%%%%
% [h2a,p2a]=ttest2(Mt2.age_days_SP_acquisition,Mt5.age_days_SP_acquisition);
% [h1b,p1b]=ttest2(Mt3.age_days_SP_acquisition,Mt5.age_days_SP_acquisition); %%%%%

%% Wilcoxon rank sum test among infants 
[pwc1,hwc1,statswc1] = ranksum(Mt1.age_days_SP_acquisition,Mt2.age_days_SP_acquisition);
[pwc2,hwc2,statswc2] = ranksum(Mt1.age_days_SP_acquisition,Mt3.age_days_SP_acquisition);
[pwc3,hwc3,statswc3] = ranksum(Mt1.age_days_SP_acquisition,Mt5.age_days_SP_acquisition);
[pwc4,hwc4,statswc4] = ranksum(Mt2.age_days_SP_acquisition,Mt3.age_days_SP_acquisition);
[pwc5,hwc5,statswc5] = ranksum(Mt2.age_days_SP_acquisition,Mt5.age_days_SP_acquisition);
[pwc6,hwc6,statswc6] = ranksum(Mt3.age_days_SP_acquisition,Mt5.age_days_SP_acquisition); 

FDR_threshold=0.05;

% % Benjamini-Hochberg (False Discovery Rate, FDR) correction 

% Store the original p-values and results in a cell array
pwc_original = [pwc1, pwc2, pwc3, pwc4, pwc5, pwc6];

adjusted_p_values_c = benjamini_hochberg_correction(pwc_original);

% Display the adjusted p-values
disp(adjusted_p_values_c);

% % Identify significant comparisons
significant_comparisons_c = find(adjusted_p_values_c <= FDR_threshold);


% Display the original p-values and their corrected values
fprintf('Original p-values\tCorrected p-values\n');
for i = 1:length(pwc_original)
    fprintf('%f\t\t\t%f\n', pwc_original(i), adjusted_p_values_c(i));
end

 
% Display significant comparisons
if ~isempty(significant_comparisons_c)
    fprintf('\nSignificant Comparisons according to Wilcoxon rank sum test (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_comparisons_c)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_comparisons_c(i), adjusted_p_values_c(significant_comparisons_c(i)));
    end
else
    disp('No significant differences were found according to Wilcoxon rank sum test with FDR correction.');
end




%%

avergae_acquisition_M1=mean(Mt1.age_days_SP_acquisition);
avergae_acquisition_M2=mean(Mt2.age_days_SP_acquisition);
avergae_acquisition_M3=mean(Mt3.age_days_SP_acquisition);
avergae_acquisition_M5=mean(Mt5.age_days_SP_acquisition);

fprintf('The average acquisition time of SP in SP -> RSV  is: %f\n',avergae_acquisition_M1 );
% Calculate the confidence interval
confidenceInterval_M1 = calculateConfidenceInterval(Mt1.age_days_SP_acquisition);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_M1(1), confidenceInterval_M1(2));


fprintf('The average acquisition time of SP in SP & RSV  is: %f\n', avergae_acquisition_M2);
% Calculate the confidence interval
confidenceInterval_M2 = calculateConfidenceInterval(Mt2.age_days_SP_acquisition);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_M2(1), confidenceInterval_M2(2));


fprintf('The average acquisition time of SP in RSV -> SP  is: %f\n', avergae_acquisition_M3);
% Calculate the confidence interval
confidenceInterval_M3 = calculateConfidenceInterval(Mt3.age_days_SP_acquisition);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_M3(1), confidenceInterval_M3(2));


fprintf('The average acquisition time of SP in SP  is: %f\n', avergae_acquisition_M5);
% Calculate the confidence interval
confidenceInterval_M5 = calculateConfidenceInterval(Mt5.age_days_SP_acquisition);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_M5(1), confidenceInterval_M5(2));



acquisition_all=[Mt1.age_days_SP_acquisition;Mt2.age_days_SP_acquisition;Mt3.age_days_SP_acquisition;Mt5.age_days_SP_acquisition];
mean_all=mean(acquisition_all);

fprintf('The average acquisition time in all infants is: %f\n', mean_all);
confidenceInterval_all = calculateConfidenceInterval(acquisition_all);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_all(1), confidenceInterval_all(2));%%





