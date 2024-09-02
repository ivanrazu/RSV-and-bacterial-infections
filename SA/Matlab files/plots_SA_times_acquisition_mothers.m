clear
clc

load('MtxGroup1_mother_with_demographic_data.mat','MtxGroup1_mother');


%%

 M=MtxGroup1_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));    
    ix=find((M.SA_Ct_Mean(indx,:))>0);    
    ix2=find((M.RSV_CT(indx,:))<99);    
    mean_SA(j)= mean( (M.SA_Ct_Mean((indx(ix)),:) ) );
    sd_SA(j)=  std( (M.SA_Ct_Mean((indx(ix)),:) ) );    
    mean_RSV(j)= mean( (M.RSV_CT((indx(ix2)),:) ) );
    %     sd_RSV(j)=  std( (M.RSV_CT((indx(ix2)),:) ) );
    mother_HIV(j)=M.Mother_HIV_status((indx(1)),:) ;
    birth_weight(j)=M.birth_weight((indx(1)),:) ;
    mother_age(j)=M.Mother_age((indx(1)),:) ;
    age_in_days_acquisition_SA(j)= time2num((M.Date_of_Visit(indx(ix(1))) - M.infant_DOB(indx(ix(1)))),"days");    
end

group=ones(1,length(subjects))*1;

Mt1=table(subjects',mean_SA',sd_SA',mean_RSV',mother_HIV',birth_weight',mother_age',group',age_in_days_acquisition_SA');
Mt1=renamevars(Mt1,"Var1","Subjects");
Mt1=renamevars(Mt1,"Var2","mean_SA");
Mt1=renamevars(Mt1,"Var3","sd_SA");
Mt1=renamevars(Mt1,"Var4","mean_RSV");
Mt1=renamevars(Mt1,"Var5","mother_HIV");
Mt1=renamevars(Mt1,"Var6","birth_weight");
Mt1=renamevars(Mt1,"Var7","mother_age");
Mt1=renamevars(Mt1,"Var8","group");
Mt1=renamevars(Mt1,"Var9","age_days_SA_acquisition");


%%

clearvars -except Mt1 Mt2

load('MtxGroup3_mother_with_demographic_data.mat','MtxGroup3_mother');


M=MtxGroup3_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));
    
     ix=find((M.SA_Ct_Mean(indx,:))>0);
     
     ix2=find((M.RSV_CT(indx,:))<99);
    
    mean_SA(j)= mean( (M.SA_Ct_Mean((indx(ix)),:) ) );
    sd_SA(j)=  std( (M.SA_Ct_Mean((indx(ix)),:) ) );
    mean_RSV(j)= mean( (M.RSV_CT((indx(ix2)),:) ) );
    mother_HIV(j)=M.Mother_HIV_status((indx(1)),:) ;
    birth_weight(j)=M.birth_weight((indx(1)),:) ;
    mother_age(j)=M.Mother_age((indx(1)),:) ;
    age_in_days_acquisition_SA(j)= time2num((M.Date_of_Visit(indx(ix(1))) - M.infant_DOB(indx(ix(1)))),"days");
        
end

group=ones(1,length(subjects))*3;

Mt3=table(subjects',mean_SA',sd_SA',mean_RSV',mother_HIV',birth_weight',mother_age',group',age_in_days_acquisition_SA');
Mt3=renamevars(Mt3,"Var1","Subjects");
Mt3=renamevars(Mt3,"Var2","mean_SA");
Mt3=renamevars(Mt3,"Var3","sd_SA");
Mt3=renamevars(Mt3,"Var4","mean_RSV");
Mt3=renamevars(Mt3,"Var5","mother_HIV");
Mt3=renamevars(Mt3,"Var6","birth_weight");
Mt3=renamevars(Mt3,"Var7","mother_age");
Mt3=renamevars(Mt3,"Var8","group");
Mt3=renamevars(Mt3,"Var9","age_days_SA_acquisition");

%%

clearvars -except Mt1 Mt2 Mt3 

load('MtxGroup4_mother_with_demographic_data.mat','MtxGroup4_mother');


M=MtxGroup4_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));
    
     ix=find((M.SA_Ct_Mean(indx,:))>0);
     
     ix2=find((M.RSV_CT(indx,:))<99);
    
    mean_SA(j)= mean( (M.SA_Ct_Mean((indx(ix)),:) ) );
    sd_SA(j)=  std( (M.SA_Ct_Mean((indx(ix)),:) ) );
    
    mean_RSV(j)= mean( (M.RSV_CT((indx(ix2)),:) ) );

    mother_HIV(j)=M.Mother_HIV_status(j,:) ;
    birth_weight(j)=M.birth_weight(j,:) ;
    
    mother_age(j)=M.Mother_age((indx(1)),:) ;
end

group=ones(1,length(subjects))*5;

Mt4=table(subjects',mean_SA',sd_SA',mean_RSV',mother_HIV',birth_weight',mother_age',group');
Mt4=renamevars(Mt4,"Var1","Subjects");
Mt4=renamevars(Mt4,"Var2","mean_SA");
Mt4=renamevars(Mt4,"Var3","sd_SA");
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
    
     ix=find((M.SA_Ct_Mean(indx,:))>0);
     
     ix2=find((M.RSV_CT(indx,:))==99);
    
    mean_SA(j)= mean( (M.SA_Ct_Mean((indx(ix)),:) ) );
    sd_SA(j)=  std( (M.SA_Ct_Mean((indx(ix)),:) ) );
    mother_HIV(j)=M.Mother_HIV_status(j,:) ;
    birth_weight(j)=M.birth_weight(j,:) ;
    mother_age(j)=M.Mother_age((indx(1)),:) ;
    age_in_days_acquisition_SA(j)= time2num((M.Date_of_Visit(indx(ix(1))) - M.infant_DOB(indx(ix(1)))),"days");

end

group=ones(1,length(subjects))*5;

mean_RSV=zeros(1,length(subjects));


Mt5=table(subjects',mean_SA',sd_SA',mean_RSV',mother_HIV',birth_weight',mother_age',group',age_in_days_acquisition_SA');
Mt5=renamevars(Mt5,"Var1","Subjects");
Mt5=renamevars(Mt5,"Var2","mean_SA");
Mt5=renamevars(Mt5,"Var3","sd_SA");
Mt5=renamevars(Mt5,"Var4","mean_RSV");
Mt5=renamevars(Mt5,"Var5","mother_HIV");
Mt5=renamevars(Mt5,"Var6","birth_weight");
Mt5=renamevars(Mt5,"Var7","mother_age");
Mt5=renamevars(Mt5,"Var8","group");
Mt5=renamevars(Mt5,"Var9","age_days_SA_acquisition");

%%

clearvars -except Mt1 Mt2 Mt3 Mt4 Mt5

load('MtxGroup6_mother_with_demographic_data.mat','MtxGroup6_mother');


M=MtxGroup6_mother;

subjects=unique(M.subject_id)';


for j=1:length(subjects)
    
    indx=find(M.subject_id==subjects(j));
    
     ix=find((M.SA_Ct_Mean(indx,:))>0);
     
     ix2=find((M.RSV_CT(indx,:))==99);
    
    mean_SA(j)= mean( (M.SA_Ct_Mean((indx(ix)),:) ) );
    sd_SA(j)=  std( (M.SA_Ct_Mean((indx(ix)),:) ) );
    mother_HIV(j)=M.Mother_HIV_status(j,:) ;
    birth_weight(j)=M.birth_weight(j,:) ;
    
    mother_age(j)=M.Mother_age((indx(1)),:) ;
end

group=ones(1,length(subjects))*6;

mean_RSV=zeros(1,length(subjects));


Mt6=table(subjects',mean_SA',sd_SA',mean_RSV',mother_HIV',birth_weight',mother_age',group');
Mt6=renamevars(Mt6,"Var1","Subjects");
Mt6=renamevars(Mt6,"Var2","mean_SA");
Mt6=renamevars(Mt6,"Var3","sd_SA");
Mt6=renamevars(Mt6,"Var4","mean_RSV");
Mt6=renamevars(Mt6,"Var5","mother_HIV");
Mt6=renamevars(Mt6,"Var6","birth_weight");
Mt6=renamevars(Mt6,"Var7","mother_age");
Mt6=renamevars(Mt6,"Var8","group");


%%
col_bl=[0 0.4470 0.7410];

X=[Mt1.age_days_SA_acquisition;Mt3.age_days_SA_acquisition;Mt5.age_days_SA_acquisition];
XX=time2num(X,"days");
Y=[ones(length(Mt1.age_days_SA_acquisition),1);2*ones(length(Mt3.age_days_SA_acquisition),1);3*ones(length(Mt5.age_days_SA_acquisition),1)];
figure
xSize = 15; Xs=xSize; ySize = 12;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

boxplot(XX,Y,'Color','k','Symbol','');

hold on
scatter(ones(length(Mt1.age_days_SA_acquisition),1),Mt1.age_days_SA_acquisition,1000,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(length(Mt3.age_days_SA_acquisition),1),Mt3.age_days_SA_acquisition,1000,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(3*ones(length(Mt5.age_days_SA_acquisition),1),Mt5.age_days_SA_acquisition,1000,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot(1,mean(Mt1.age_days_SA_acquisition),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
plot(2,mean(Mt3.age_days_SA_acquisition),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
plot(3,mean(Mt5.age_days_SA_acquisition),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
xlim([0.5,3.5])
ylim([-5,120])
xticks((1:1:3))
set(gca,'Fontsize',60);box on;
set(gca, 'XTickLabel', {'SA$\rightarrow$RSV','RSV$\rightarrow$SA','SA'},'TickLabelInterpreter','latex','Fontsize',43)
ylabel('SA acquisition (days)','Interpreter','latex')
set(gca,'Fontsize',60);box on;
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(30)
ylim([-5,125])
set(gca,'Fontsize',60);box on;

%% test for normality

[hc1,pc1] = kstest(Mt1.age_days_SA_acquisition);
[hc3,pc3] = kstest(Mt3.age_days_SA_acquisition);
[hc4,pc4] = kstest(Mt5.age_days_SA_acquisition);

pc=[pc1,pc3,pc4];

% Check for significant differences
alpha = 0.05;  % Set significance level
for i=1:3
    if pc(i) < alpha
        fprintf('Mothers Group %d is NOT normally distributed.\n', i );
    else
        fprintf('Mothers Group %d is normally distributed.\n', i );
    end
end


%% T-test NOT recommended as samples are NOT normally distributed  
% [h2,p2]=ttest2(Mt1.age_days_SA_acquisition,Mt3.age_days_SA_acquisition); %%%%%%%  p=0.0377
% [h3,p3]=ttest2(Mt1.age_days_SA_acquisition,Mt5.age_days_SA_acquisition) ;
% [h1b,p1b]=ttest2(Mt3.age_days_SA_acquisition,Mt5.age_days_SA_acquisition) ;


%% Wilcoxon rank sum test among infants 
[pwc2,hwc2,statswc2] = ranksum(Mt1.age_days_SA_acquisition,Mt3.age_days_SA_acquisition); 
[pwc3,hwc3,statswc3] = ranksum(Mt1.age_days_SA_acquisition,Mt5.age_days_SA_acquisition); 
[pwc6,hwc6,statswc6] = ranksum(Mt3.age_days_SA_acquisition,Mt5.age_days_SA_acquisition); 

FDR_threshold=0.05;

% % Benjamini-Hochberg (False Discovery Rate, FDR) correction 

% Store the original p-values and results in a cell array
pwc_original = [pwc2, pwc3, pwc6];

adjusted_p_values_m = benjamini_hochberg_correction(pwc_original);

% Display the adjusted p-values
disp(adjusted_p_values_m);

% % Identify significant comparisons
significant_comparisons_c = find(adjusted_p_values_m <= FDR_threshold);


% Display the original p-values and their corrected values
fprintf('Original p-values\tCorrected p-values\n');
for i = 1:length(pwc_original)
    fprintf('%f\t\t\t%f\n', pwc_original(i), adjusted_p_values_m(i));
end

 
% Display significant comparisons
if ~isempty(significant_comparisons_c)
    fprintf('\nSignificant Comparisons according to Wilcoxon rank sum test (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_comparisons_c)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_comparisons_c(i), adjusted_p_values_m(significant_comparisons_c(i)));
    end
else
    disp('No significant differences were found according to Wilcoxon rank sum test with FDR correction.');
end


%%

avergae_acquisition_M1=mean(Mt1.age_days_SA_acquisition);
avergae_acquisition_M3=mean(Mt3.age_days_SA_acquisition);
avergae_acquisition_M5=mean(Mt5.age_days_SA_acquisition);

avergae_acquisition_M1_3=mean([Mt1.age_days_SA_acquisition;Mt3.age_days_SA_acquisition]);


fprintf('The average acquisition time of SA in SA -> RSV  is: %f\n',avergae_acquisition_M1 );
% Calculate the confidence interval
confidenceInterval_M1 = calculateConfidenceInterval(Mt1.age_days_SA_acquisition);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_M1(1), confidenceInterval_M1(2));

fprintf('The average acquisition time of SA in RSV -> SA  is: %f\n', avergae_acquisition_M3);
% Calculate the confidence interval
confidenceInterval_M3 = calculateConfidenceInterval(Mt3.age_days_SA_acquisition);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_M3(1), confidenceInterval_M3(2));


fprintf('The average acquisition time of SA in SA  is: %f\n', avergae_acquisition_M5);
% Calculate the confidence interval
confidenceInterval_M5 = calculateConfidenceInterval(Mt5.age_days_SA_acquisition);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_M5(1), confidenceInterval_M5(2));


acquisition_all=[Mt1.age_days_SA_acquisition;Mt3.age_days_SA_acquisition;Mt5.age_days_SA_acquisition];
mean_all=mean(acquisition_all);

fprintf('The average acquisition time in all infants is: %f\n', mean_all);
confidenceInterval_all = calculateConfidenceInterval(acquisition_all);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_all(1), confidenceInterval_all(2));%%



fprintf('The average SA acquisition time for infants in SA->RSV,SA&RSV,RSV->SA is: %f\n',avergae_acquisition_M1_3)
% Calculate the confidence interval
confidenceInterval_M1_3 = calculateConfidenceInterval([Mt1.age_days_SA_acquisition;Mt3.age_days_SA_acquisition]);
fprintf('%.2f%% Confidence Interval: [%.2f, %.2f]\n', 95, confidenceInterval_M1_3(1), confidenceInterval_M1_3(2));





