%Plots boxplots of SA CT on Mothers only, infants only, and both combined
clear
clc
warning off

dot_size=400;

load('MtxGroup1_child_with_demographic_data.mat','MtxGroup1_child');
load('MtxGroup2_child_with_demographic_data.mat','MtxGroup2_child');
load('MtxGroup3_child_with_demographic_data.mat','MtxGroup3_child');
load('MtxGroup4_child_with_demographic_data.mat','MtxGroup4_child');
load('MtxGroup5_child_with_demographic_data.mat','MtxGroup5_child');


load('MtxGroup1_mother_with_demographic_data.mat','MtxGroup1_mother');
load('MtxGroup3_mother_with_demographic_data.mat','MtxGroup3_mother');
load('MtxGroup4_mother_with_demographic_data.mat','MtxGroup4_mother');
load('MtxGroup5_mother_with_demographic_data.mat','MtxGroup5_mother');


M1m=MtxGroup1_mother;
M3m=MtxGroup3_mother;
M4m=MtxGroup4_mother;
M5m=MtxGroup5_mother;


M1c=MtxGroup1_child;
M2c=MtxGroup2_child;
M3c=MtxGroup3_child;
M4c=MtxGroup4_child;
M5c=MtxGroup5_child;
%%

SA_positive_only_M1c=M1c.SA_Ct_Mean(M1c.SA_Ct_Mean>0);
mean_SA_positive_M1c=mean(SA_positive_only_M1c);

SA_positive_only_M2c=M2c.SA_Ct_Mean(M2c.SA_Ct_Mean>0);
mean_SA_positive_M2c=mean(SA_positive_only_M2c);

SA_positive_only_M3c=M3c.SA_Ct_Mean(M3c.SA_Ct_Mean>0);
mean_SA_positive_M3c=mean(SA_positive_only_M3c);

SA_positive_only_M5c=M5c.SA_Ct_Mean(M5c.SA_Ct_Mean>0);
mean_SA_positive_M5c=mean(SA_positive_only_M5c);


SA_positive_only_M1m=M1m.SA_Ct_Mean(M1m.SA_Ct_Mean>0);
mean_SA_positive_M1m=mean(SA_positive_only_M1m);

SA_positive_only_M3m=M3m.SA_Ct_Mean(M3m.SA_Ct_Mean>0);
mean_SA_positive_M3m=mean(SA_positive_only_M3m);

SA_positive_only_M5m=M5m.SA_Ct_Mean(M5m.SA_Ct_Mean>0);
mean_SA_positive_M5m=mean(SA_positive_only_M5m);
    

MM1c=table(SA_positive_only_M1c, ones(1,length(SA_positive_only_M1c))');
MM2c=table(SA_positive_only_M2c,2*ones(1,length(SA_positive_only_M2c))');
MM3c=table(SA_positive_only_M3c,3*ones(1,length(SA_positive_only_M3c))');
MM5c=table(SA_positive_only_M5c,4*ones(1,length(SA_positive_only_M5c))');

MM1c=renamevars(MM1c,["Var2","SA_positive_only_M1c"],["group","SA_Ct"]);
MM2c=renamevars(MM2c,["Var2","SA_positive_only_M2c"],["group","SA_Ct"]);
MM3c=renamevars(MM3c,["Var2","SA_positive_only_M3c"],["group","SA_Ct"]);
MM5c=renamevars(MM5c,["Var2","SA_positive_only_M5c"],["group","SA_Ct"]);


MM1m=table(SA_positive_only_M1m, ones(1,length(SA_positive_only_M1m))');
MM3m=table(SA_positive_only_M3m,3*ones(1,length(SA_positive_only_M3m))');
MM5m=table(SA_positive_only_M5m,4*ones(1,length(SA_positive_only_M5m))');

MM1m=renamevars(MM1m,["Var2","SA_positive_only_M1m"],["group","SA_Ct"]);
MM3m=renamevars(MM3m,["Var2","SA_positive_only_M3m"],["group","SA_Ct"]);
MM5m=renamevars(MM5m,["Var2","SA_positive_only_M5m"],["group","SA_Ct"]);


MMm=[MM1m;MM3m;MM5m];
meanMMm = groupsummary(MMm,'group','mean');

MMc=[MM1c;MM3c;MM5c];
meanMMc = groupsummary(MMc,'group','mean');


%% boxplot mothers vs infants using boxplotGroup

n1=max([length(MM1m.SA_Ct),length(MM3m.SA_Ct),length(MM5m.SA_Ct)]);
n2=max([length(MM1c.SA_Ct),length(MM3c.SA_Ct),length(MM5c.SA_Ct)]);


MM1m.SA_Ct(end+1:n1)=nan;
MM3m.SA_Ct(end+1:n1)=nan;
MM5m.SA_Ct(end+1:n1)=nan;

MM1c.SA_Ct(end+1:n2)=nan;
MM2c.SA_Ct(end+1:n2)=nan;
MM3c.SA_Ct(end+1:n2)=nan;
MM5c.SA_Ct(end+1:n2)=nan;

X={[MM1m.SA_Ct,MM3m.SA_Ct,MM5m.SA_Ct],[MM1c.SA_Ct,MM3c.SA_Ct,MM5c.SA_Ct]};

col_bl=[0 0.4470 0.7410];
col_blk=[0,0,0];

%%
cols=1;
rows=1;
fig=figure;
xSize = cols*18; Xs=xSize; ySize = rows*12;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[100 25 xSize*50 ySize*55]);

boxplotGroup(X,'Colors',[col_bl;col_blk],'GroupType','betweenGroups','Symbol','','groupLines', true);

hold on;
scatter(1*ones(1,length(MM1m.SA_Ct)),MM1m.SA_Ct,200,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(1,length(MM1c.SA_Ct)),MM1c.SA_Ct,200,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(4*ones(1,length(MM3m.SA_Ct)),MM3m.SA_Ct,200,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(5*ones(1,length(MM3c.SA_Ct)),MM3c.SA_Ct,200,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(7*ones(1,length(MM5m.SA_Ct)),MM5m.SA_Ct,200,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(8*ones(1,length(MM5c.SA_Ct)),MM5c.SA_Ct,200,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot([1,4,7],meanMMm.mean_SA_Ct,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
hold on
plot([2,5,8],meanMMc.mean_SA_Ct,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
set(gca,'Fontsize',60);box on;
xticks([1.5,4.5,7.5])
ylim([18,45])
set(gca, 'XTickLabel', {'SA$\rightarrow$RSV','RSV$\rightarrow$SA','SA'},'TickLabelInterpreter','latex','Fontsize',60)
ylabel('SA Ct','Interpreter','latex','Fontsize',60)
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(20)
xline(3,'LineWidth',3)
xline(6,'LineWidth',3)


% % plots line to denote significant difference
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.1)])
% xt = get(gca, 'XTick');
xt=[1 2];

yxis_reverse=1;

if yxis_reverse==0
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.05, 'k-','LineWidth',6)
ylim([18,45])
end

if yxis_reverse==1
    set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*0.45, 'k-','LineWidth',6)
end

dim1 = [.245 .62 .9 .3];

name={'$ * $'};
annotation('textbox',dim1,'String',name,'interpreter','latex','Fontsize',55,'Color', 'k','EdgeColor','none');

set(gca,'Fontsize',60);box on;


%% two sample Wilcoxon rank sum test Infants vs Mothers by group

[p1a,h1a]=ranksum(MM1m.SA_Ct,MM1c.SA_Ct);%%% p=0.010464.
[p3a,h3a]=ranksum(MM3m.SA_Ct,MM3c.SA_Ct) ;
[p4a,h4a]=ranksum(MM5m.SA_Ct,MM5c.SA_Ct) ;

pmc=[p1a,p3a,p4a];

% Check for significant differences
alpha = 0.05;  % Set significance level
group_labels =  {'SA->RSV', 'RSV->SA', 'SA only'};

fprintf('SA&RSV was not considered, not enough data.\n');

for i=1:3
    if pmc(i) < alpha
        fprintf('Infants-Mothers are significantly different in %s p=%f.\n', group_labels{i}, pmc(i) );
    else
        fprintf('Infants-Mothers are NOT significantly different in %s p=%f.\n', group_labels{i}, pmc(i) );
    end
end


%% boxplot mothers only

XXm=[MM1m.SA_Ct;MM3m.SA_Ct;MM5m.SA_Ct];
YYm=[ones(size(MM1m.SA_Ct,1),1);2*ones(size(MM1m.SA_Ct,1),1);3*ones(size(MM1m.SA_Ct,1),1)];

figure;
xSize = 15; Xs=xSize; ySize = 13.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

boxplot(XXm,YYm,'Color','k','Symbol','');

hold on;
scatter(1*ones(1,length(MM1m.SA_Ct)),MM1m.SA_Ct,dot_size,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(1,length(MM3m.SA_Ct)),MM3m.SA_Ct,dot_size,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(3*ones(1,length(MM5m.SA_Ct)),MM5m.SA_Ct,dot_size,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot([1,2,3],meanMMm.mean_SA_Ct,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',15);
set(gca,'Fontsize',60);box on;
xticks([1,2,3])
ylim([18,45])
set(gca, 'XTickLabel', {'SA$\rightarrow$RSV','RSV$\rightarrow$SA','SA'},'TickLabelInterpreter','latex','Fontsize',50)
ylabel('SA Ct','Interpreter','latex','Fontsize',60)
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(30)


% % plots line to denote significant difference
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.1)])
xt=[1,2,3];

if yxis_reverse==0
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.02, 'k-','LineWidth',6)
    hold on
    plot(xt([2 3]), [1 1]*max(yt)*1.05, 'k-','LineWidth',6)
    hold on
    plot(xt([1 3]), [1 1]*max(yt)*1.09, 'k-','LineWidth',6)
    dim1 = [.37 .48 .9 .3];
    dim2 = [.63 .50 .9 .3];
    dim3 = [.5 .55 .9 .3];
end

set(gca,'Fontsize',60);box on;

ylim([18,50])

if yxis_reverse==1
    set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*0.53, 'k-','LineWidth',6)
    hold on
    plot(xt([2 3]), [1 1]*max(yt)*0.44, 'k-','LineWidth',6)
    hold on
    plot(xt([1 3]), [1 1]*max(yt)*0.49, 'k-','LineWidth',6)
    hold on
    dim1 = [.37 .53 .9 .3];
    dim2 = [.63 .62 .9 .3];
    dim3 = [.5 .57 .9 .3];

    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
end


name={'$ * $'};
annotation('textbox',dim1,'String',name,'interpreter','latex','Fontsize',55,'Color', 'k','EdgeColor','none');
hold on
annotation('textbox',dim2,'String',name,'interpreter','latex','Fontsize',55,'Color', 'k','EdgeColor','none');
hold on
annotation('textbox',dim3,'String',name,'interpreter','latex','Fontsize',55,'Color', 'k','EdgeColor','none');

%%  t test among mother groups

[h1m,p1m]=ttest2(MM1m.SA_Ct,MM3m.SA_Ct);
[h2m,p2m]=ttest2(MM1m.SA_Ct,MM5m.SA_Ct);
[h3m,p3m]=ttest2(MM3m.SA_Ct,MM5m.SA_Ct);


pim=[p1m,p2m,p3m];

% Find which p-values are lower than alpha
significant_indices_ttest_m = find(pim < alpha);

% Create a cell array of variable names
variable_names_m = {'p1m', 'p2m', 'p3m'};

% Display the actual p-values and their variable names
if ~isempty(significant_indices_ttest_m)
    disp('There are significant differences among the mother groups according to t-test.');
    disp('Significant p-values:');
    for i = significant_indices_ttest_m
        variable_name_m = variable_names_m{i};
        p_value_m = pim(i);
        disp([variable_name_m, ' = ', num2str(p_value_m)]);
    end
else
    disp('No significant differences among the mother groups were found according to t-test.');
end
 
%% % Benjamini-Hochberg (False Discovery Rate, FDR) correction for t-test

FDR_threshold=0.05;

adjusted_p_values_tstest_Mothers = benjamini_hochberg_correction(pim);

% % Identify significant comparisons
significant_indices_ttest_BH_m = find(adjusted_p_values_tstest_Mothers <= FDR_threshold);

group_comparisons={'SA->RSV	vs RSV->SA ' , 'SA->RSV vs SA   ',   'RSV->SA vs SA   '};


% Display the original p-values and their corrected values
fprintf('T-test Mothers \n');
fprintf('Comparison\tOriginal p-values\tBH-Corrected p-values\n');
for i = 1:length(pim)
    fprintf('%s\t%f\t\t\t%f\n',group_comparisons{i} ,pim(i), adjusted_p_values_tstest_Mothers(i));
end


% Display significant comparisons
if ~isempty(significant_indices_ttest_BH_m)
    fprintf('\nSignificant Comparisons in Mother Groups according to t-test (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_indices_ttest_BH_m)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_indices_ttest_BH_m(i), adjusted_p_values_tstest_Mothers(significant_indices_ttest_BH_m(i)));
    end
else
    disp('No significant differences in Mother Groups according to t-test were found after FDR correction.');
end

%% Wilcoxon rank sum test among momthers 
[pwm1,hwm1,statswm1] = ranksum(MM1m.SA_Ct,MM3m.SA_Ct);
[pwm2,hwm2,statswm2] = ranksum(MM1m.SA_Ct,MM5m.SA_Ct); 
[pwm3,hwm3,statswm3] = ranksum(MM3m.SA_Ct,MM5m.SA_Ct);


FDR_threshold=0.05;

% % Benjamini-Hochberg (False Discovery Rate, FDR) correction 

% p-values obtained from Wilcoxon rank-sum tests
pwm_original = [pwm1, pwm2, pwm3];


adjusted_p_values_WRS_Mothers = benjamini_hochberg_correction(pwm_original);

% % Identify significant comparisons
significant_indices_WRS_test_BH_m = find(adjusted_p_values_WRS_Mothers <= FDR_threshold);


% Display the original p-values and their corrected values
fprintf('WRS-test Mothers \n');
fprintf('Comparison\tOriginal p-values\tBH-Corrected p-values\n');
for i = 1:length(pwm_original)
    fprintf('%s\t%f\t\t\t%f\n',group_comparisons{i}, pwm_original(i), adjusted_p_values_WRS_Mothers(i));
end

 
% Display significant comparisons
if ~isempty(significant_indices_WRS_test_BH_m)
    fprintf('\nSignificant Comparisons (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_indices_WRS_test_BH_m)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_indices_WRS_test_BH_m(i), adjusted_p_values_WRS_Mothers(significant_indices_WRS_test_BH_m(i)));
    end
else
    disp('No significant differences were found after FDR correction.');
end



%% boxplot infants only

MMc=[MM1c;MM2c;MM3c;MM5c];
meanMMc = groupsummary(MMc,'group','mean');

XXc=[MM1c.SA_Ct;MM2c.SA_Ct;MM3c.SA_Ct;MM5c.SA_Ct];
YYc=[ones(size(MM1c.SA_Ct,1),1);2*ones(size(MM1c.SA_Ct,1),1);3*ones(size(MM1c.SA_Ct,1),1);4*ones(size(MM1c.SA_Ct,1),1)];


cols=1;
rows=1;
figure
xSize = 15; Xs=xSize; ySize = 13.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);
boxplot(XXc,YYc,'Color','k','Symbol','');
set(gca,'Fontsize',60);box on;
hold on;
scatter(1*ones(1,length(MM1c.SA_Ct)),MM1c.SA_Ct,dot_size,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(1,length(MM2c.SA_Ct)),MM2c.SA_Ct,dot_size,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(3*ones(1,length(MM3c.SA_Ct)),MM3c.SA_Ct,dot_size,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(4*ones(1,length(MM5c.SA_Ct)),MM5c.SA_Ct,dot_size,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot([0,1,2,3,4],meanMMc.mean_SA_Ct,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20);
set(gca,'Fontsize',40);box on;
xticks([1,2,3,4])
ylim([18,45])
set(gca, 'XTickLabel', {'SA$\rightarrow$RSV','SA\&RSV','RSV$\rightarrow$SA','SA'},'TickLabelInterpreter','latex','FontSize',25)
ylabel('SA Ct','Interpreter','latex')
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(30)
ylim([18,50])
set(gca,'Fontsize',60);box on;


if yxis_reverse==1
    set(gca, 'YDir', 'reverse')
    ylim([15,45])
    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
end

%% t-test among infant groups
[h1c,p1c]=ttest2(MM1c.SA_Ct,MM2c.SA_Ct);
[h2c,p2c]=ttest2(MM1c.SA_Ct,MM3c.SA_Ct);
[h3c,p3c]=ttest2(MM1c.SA_Ct,MM5c.SA_Ct);
[h4c,p4c]=ttest2(MM2c.SA_Ct,MM3c.SA_Ct);
[h5c,p5c]=ttest2(MM2c.SA_Ct,MM5c.SA_Ct);
[h6c,p6c]=ttest2(MM3c.SA_Ct,MM5c.SA_Ct);

pic=[p1c,p2c,p3c,p4c,p5c,p6c];

% Find which p-values are lower than alpha
significant_indices_ttest_c = find(pic < alpha);

% Create a cell array of variable names
variable_names_c = {'p1c', 'p2c', 'p3c', 'p4c', 'p5c', 'p6c'};

% Display the actual p-values and their variable names
if ~isempty(significant_indices_ttest_c)
    disp('There are significant differences among the infant groups according to t-test.');
    disp('Significant p-values:');
    for i = significant_indices_ttest_c
        variable_name_c = variable_names_c{i};
        p_value_c = pic(i);
        disp([variable_name_c, ' = ', num2str(p_value_c)]);
    end
else
    disp('No significant differences among the infant groups were found according to t-test.');
end


%% Benjamini-Hochberg (False Discovery Rate, FDR) correction 

FDR_threshold=0.05;

% Store the original p-values and results in a cell array
pct_original = [p1c, p2c, p3c, p4c, p5c, p6c];

adjusted_p_values_ttest_Infants = benjamini_hochberg_correction(pct_original);


% % Identify significant comparisons
significant_indices_ttest_BH_c = find(adjusted_p_values_ttest_Infants <= FDR_threshold);

group_comparisons={'SA->RSV	vs SA&RSV  ', 'SA->RSV	vs RSV->SA ' , 'SA->RSV vs SA   ', 'SA&RSV  vs RSV->SA ', 'SA&RSV vs SA     ', 'RSV->SA vs SA   '};


% Display the original p-values and their corrected values
fprintf('T-test Infants\n');
fprintf('Comparison\tOriginal p-values\tCorrected p-values\n');
for i = 1:length(pct_original)
    fprintf('%s\t%f\t\t\t%f\n',group_comparisons{i}, pct_original(i), adjusted_p_values_ttest_Infants(i));
end


% Display significant comparisons
if ~isempty(significant_indices_ttest_BH_c)
    fprintf('\nSignificant Comparisons in Infant Groups according to t-test (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_indices_ttest_BH_c)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_indices_ttest_BH_c(i), adjusted_p_values_ttest_Infants(significant_indices_ttest_BH_c(i)));
    end
else
    disp('No significant differences in Infant Groups were found according to t-test with FDR correction.');
end


%% Wilcoxon rank sum test among infants 
[pwc1,hwc1,statswc1] = ranksum(MM1c.SA_Ct,MM2c.SA_Ct);
[pwc2,hwc2,statswc2] = ranksum(MM1c.SA_Ct,MM3c.SA_Ct); 
[pwc3,hwc3,statswc3] = ranksum(MM1c.SA_Ct,MM5c.SA_Ct);
[pwc4,hwc4,statswc4] = ranksum(MM2c.SA_Ct,MM3c.SA_Ct);
[pwc5,hwc5,statswc5] = ranksum(MM2c.SA_Ct,MM5c.SA_Ct);
[pwc6,hwc6,statswc6] = ranksum(MM3c.SA_Ct,MM5c.SA_Ct); 

FDR_threshold=0.05;

% % Benjamini-Hochberg (False Discovery Rate, FDR) correction 

% Store the original p-values and results in a cell array
pwc_original = [pwc1, pwc2, pwc3, pwc4, pwc5, pwc6];

adjusted_p_values_WRS_Infants = benjamini_hochberg_correction(pwc_original);


% % Identify significant comparisons
significant_indices_WRS_BH_test_c = find(adjusted_p_values_WRS_Infants <= FDR_threshold);


% Display the original p-values and their corrected values
fprintf('WRS-test Infants\n');
fprintf('Comparison\tOriginal p-values\tCorrected p-values\n');
for i = 1:length(pwc_original)
    fprintf('%s\t%f\t\t\t%f\n',group_comparisons{i}, pwc_original(i), adjusted_p_values_WRS_Infants(i));
end

 
% Display significant comparisons
if ~isempty(significant_indices_WRS_BH_test_c)
    fprintf('\nSignificant Comparisons (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_indices_WRS_BH_test_c)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_indices_WRS_BH_test_c(i), adjusted_p_values_WRS_Infants(significant_indices_WRS_BH_test_c(i)));
    end
else
    disp('No significant differences were found after FDR correction.');
end

