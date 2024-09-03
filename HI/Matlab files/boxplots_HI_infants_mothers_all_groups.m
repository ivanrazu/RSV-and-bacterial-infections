%Plots boxplots of HI CT on Mothers only, infants only, and both combined
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
load('MtxGroup2_mother_with_demographic_data.mat','MtxGroup2_mother');
load('MtxGroup3_mother_with_demographic_data.mat','MtxGroup3_mother');
load('MtxGroup4_mother_with_demographic_data.mat','MtxGroup4_mother');
load('MtxGroup5_mother_with_demographic_data.mat','MtxGroup5_mother');

M1m=MtxGroup1_mother;
M2m=MtxGroup2_mother;
M3m=MtxGroup3_mother;
M4m=MtxGroup4_mother;
M5m=MtxGroup5_mother;

M1c=MtxGroup1_child;
M2c=MtxGroup2_child;
M3c=MtxGroup3_child;
M4c=MtxGroup4_child;
M5c=MtxGroup5_child;
%%

HI_positive_only_M1c=M1c.HI_Ct_Mean(M1c.HI_Ct_Mean>0);
mean_HI_positive_M1c=mean(HI_positive_only_M1c);

HI_positive_only_M2c=M2c.HI_Ct_Mean(M2c.HI_Ct_Mean>0);
mean_HI_positive_M2c=mean(HI_positive_only_M2c);

HI_positive_only_M3c=M3c.HI_Ct_Mean(M3c.HI_Ct_Mean>0);
mean_HI_positive_M3c=mean(HI_positive_only_M3c);

HI_positive_only_M5c=M5c.HI_Ct_Mean(M5c.HI_Ct_Mean>0);
mean_HI_positive_M5c=mean(HI_positive_only_M5c);


HI_positive_only_M1m=M1m.HI_Ct_Mean(M1m.HI_Ct_Mean>0);
mean_HI_positive_M1m=mean(HI_positive_only_M1m);

HI_positive_only_M2m=M2m.HI_Ct_Mean(M2m.HI_Ct_Mean>0);
mean_HI_positive_M2m=mean(HI_positive_only_M2m);

HI_positive_only_M3m=M3m.HI_Ct_Mean(M3m.HI_Ct_Mean>0);
mean_HI_positive_M3m=mean(HI_positive_only_M3m);

HI_positive_only_M5m=M5m.HI_Ct_Mean(M5m.HI_Ct_Mean>0);
mean_HI_positive_M5m=mean(HI_positive_only_M5m);
    

MM1c=table(HI_positive_only_M1c, ones(1,length(HI_positive_only_M1c))');
MM2c=table(HI_positive_only_M2c,2*ones(1,length(HI_positive_only_M2c))');
MM3c=table(HI_positive_only_M3c,3*ones(1,length(HI_positive_only_M3c))');
MM5c=table(HI_positive_only_M5c,4*ones(1,length(HI_positive_only_M5c))');

MM1c=renamevars(MM1c,["Var2","HI_positive_only_M1c"],["group","HI_Ct"]);
MM2c=renamevars(MM2c,["Var2","HI_positive_only_M2c"],["group","HI_Ct"]);
MM3c=renamevars(MM3c,["Var2","HI_positive_only_M3c"],["group","HI_Ct"]);
MM5c=renamevars(MM5c,["Var2","HI_positive_only_M5c"],["group","HI_Ct"]);


MM1m=table(HI_positive_only_M1m, ones(1,length(HI_positive_only_M1m))');
MM2m=table(HI_positive_only_M2m,2*ones(1,length(HI_positive_only_M2m))');
MM3m=table(HI_positive_only_M3m,3*ones(1,length(HI_positive_only_M3m))');
MM5m=table(HI_positive_only_M5m,4*ones(1,length(HI_positive_only_M5m))');

MM1m=renamevars(MM1m,["Var2","HI_positive_only_M1m"],["group","HI_Ct"]);
MM2m=renamevars(MM2m,["Var2","HI_positive_only_M2m"],["group","HI_Ct"]);
MM3m=renamevars(MM3m,["Var2","HI_positive_only_M3m"],["group","HI_Ct"]);
MM5m=renamevars(MM5m,["Var2","HI_positive_only_M5m"],["group","HI_Ct"]);

MMm=[MM1m;MM2m;MM3m;MM5m];
meanMMm = groupsummary(MMm,'group','mean');

MHI=[MM1c;MM2c;MM3c;MM5c];
meanMHI = groupsummary(MHI,'group','mean');

%% boxplot mothers vs infants using boxplotGroup

x={MM1m.HI_Ct,MM2m.HI_Ct,MM3m.HI_Ct,MM5m.HI_Ct};

n1=max([length(MM1m.HI_Ct),length(MM2m.HI_Ct),length(MM3m.HI_Ct),length(MM5m.HI_Ct)]);
n2=max([length(MM1c.HI_Ct),length(MM2c.HI_Ct),length(MM3c.HI_Ct),length(MM5c.HI_Ct)]);


MM1m.HI_Ct(end+1:n1)=nan;
MM2m.HI_Ct(end+1:n1)=nan;
MM3m.HI_Ct(end+1:n1)=nan;
MM5m.HI_Ct(end+1:n1)=nan;

MM1c.HI_Ct(end+1:n2)=nan;
MM2c.HI_Ct(end+1:n2)=nan;
MM3c.HI_Ct(end+1:n2)=nan;
MM5c.HI_Ct(end+1:n2)=nan;

X={[MM1m.HI_Ct,MM2m.HI_Ct,MM3m.HI_Ct,MM5m.HI_Ct],[MM1c.HI_Ct,MM2c.HI_Ct,MM3c.HI_Ct,MM5c.HI_Ct]};

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
scatter(1*ones(1,length(MM1m.HI_Ct)),MM1m.HI_Ct,200,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(1,length(MM1c.HI_Ct)),MM1c.HI_Ct,200,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(4*ones(1,length(MM2m.HI_Ct)),MM2m.HI_Ct,200,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(5*ones(1,length(MM2c.HI_Ct)),MM2c.HI_Ct,200,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(7*ones(1,length(MM3m.HI_Ct)),MM3m.HI_Ct,200,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(8*ones(1,length(MM3c.HI_Ct)),MM3c.HI_Ct,200,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(10*ones(1,length(MM5m.HI_Ct)),MM5m.HI_Ct,200,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(11*ones(1,length(MM5c.HI_Ct)),MM5c.HI_Ct,200,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot([1,4,7,10],meanMMm.mean_HI_Ct,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
hold on
plot([2,5,8,11],meanMHI.mean_HI_Ct,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10);
set(gca,'Fontsize',60);box on;
xticks([1.5,4.5,7.5,10.5])
ylim([18,45])
set(gca, 'XTickLabel', {'HI$\rightarrow$RSV','HI\&RSV','RSV$\rightarrow$HI','HI'},'TickLabelInterpreter','latex','Fontsize',60)
ylabel('HI Ct','Interpreter','latex','FontSize',60)
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(20)
xline(3,'LineWidth',3)
xline(6,'LineWidth',3)
xline(6,'LineWidth',3)
xline(9,'LineWidth',3)

set(gca,'Fontsize',60);box on;
yxis_reverse=1;

if yxis_reverse==1
set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
end

%% %% two sample Wilcoxon rank sum test Infants vs Mothers by group
 % t-test
% [h1a,p1a]=ttest2(MM1m.HI_Ct,MM1c.HI_Ct);
% [h2a,p2a]=ttest2(MM2m.HI_Ct,MM2c.HI_Ct);
% [h3a,p3a]=ttest2(MM3m.HI_Ct,MM3c.HI_Ct);
% [h4a,p4a]=ttest2(MM5m.HI_Ct,MM5c.HI_Ct);

[p1a,h1a]=ranksum(MM1m.HI_Ct,MM1c.HI_Ct);
[p2a,h2a]=ranksum(MM2m.HI_Ct,MM2c.HI_Ct);
[p3a,h3a]=ranksum(MM3m.HI_Ct,MM3c.HI_Ct);
[p4a,h4a]=ranksum(MM5m.HI_Ct,MM5c.HI_Ct);
%%
pmc=[p1a,p2a,p3a,p4a];

% Check for significant differences
alpha = 0.05;  % Set significance level
for i=1:4
    if pmc(i) < alpha
        fprintf('Infants-Mothers Group %d are significantly different p=%f.\n', i, pmc(i) );
    else
        fprintf('Infants-Mothers Group %d are NOT significantly different p=%f.\n', i, pmc(i) );
    end
end

%% boxplot mothers only


XXm=[MM1m.HI_Ct;MM2m.HI_Ct;MM3m.HI_Ct;MM5m.HI_Ct];
YYm=[ones(size(MM1m.HI_Ct,1),1);2*ones(size(MM1m.HI_Ct,1),1);3*ones(size(MM1m.HI_Ct,1),1);4*ones(size(MM1m.HI_Ct,1),1)];

figure;
xSize = 15; Xs=xSize; ySize = 13.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

boxplot(XXm,YYm,'Color','k','Symbol','');
set(gca,'Fontsize',60);box on;

hold on;
scatter(1*ones(1,length(MM1m.HI_Ct)),MM1m.HI_Ct,dot_size,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(1,length(MM2m.HI_Ct)),MM2m.HI_Ct,dot_size,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(3*ones(1,length(MM3m.HI_Ct)),MM3m.HI_Ct,dot_size,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(4*ones(1,length(MM5m.HI_Ct)),MM5m.HI_Ct,dot_size,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot([1,2,3,4],meanMMm.mean_HI_Ct,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20);
set(gca,'Fontsize',40);box on;
xticks([1,2,3,4])
ylim([18,50])
set(gca, 'XTickLabel', {'HI$\rightarrow$RSV','HI\&RSV','RSV$\rightarrow$HI','HI'},'TickLabelInterpreter','latex','Fontsize',50)
ylabel('HI Ct','Interpreter','latex')
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(30)

set(gca,'Fontsize',60);box on;
ylim([18,50])

if yxis_reverse==1
    set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
end

%%  t test among mother groups

[hm1,p1m]=ttest2(MM1m.HI_Ct,MM2m.HI_Ct);
[hm2,p2m]=ttest2(MM1m.HI_Ct,MM3m.HI_Ct);
[hm3,p3m]=ttest2(MM1m.HI_Ct,MM5m.HI_Ct);  %%%%%  0.0113
[hm4,p4m]=ttest2(MM2m.HI_Ct,MM3m.HI_Ct);
[hm5,p5m]=ttest2(MM2m.HI_Ct,MM5m.HI_Ct);
[hm6,p6m]=ttest2(MM3m.HI_Ct,MM5m.HI_Ct);

pim=[p1m,p2m,p3m,p4m,p5m,p6m];

% Find which p-values are lower than alpha
significant_indices_ttest_m = find(pim < alpha);

% Create a cell array of variable names
variable_names_m = {'p1m', 'p2m', 'p3m', 'p4m', 'p5m', 'p6m'};

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

group_comparisons={'HI->RSV	vs HI&RSV  ', 'HI->RSV	vs RSV->HI ' , 'HI->RSV vs HI   ', 'HI&RSV  vs RSV->HI ', 'HI&RSV vs HI     ', 'RSV->HI vs HI   '};


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
[pwm1,hwm1,statswm1] = ranksum(MM1m.HI_Ct,MM2m.HI_Ct);
[pwm2,hwm2,statswm2] = ranksum(MM1m.HI_Ct,MM3m.HI_Ct); 
[pwm3,hwm3,statswm3] = ranksum(MM1m.HI_Ct,MM5m.HI_Ct);
[pwm4,hwm4,statswm4] = ranksum(MM2m.HI_Ct,MM3m.HI_Ct);
[pwm5,hwm5,statswm5] = ranksum(MM2m.HI_Ct,MM5m.HI_Ct);
[pwm6,hwm6,statswm6] = ranksum(MM3m.HI_Ct,MM5m.HI_Ct); 

FDR_threshold=0.05;

% % Benjamini-Hochberg (False Discovery Rate, FDR) correction 

% p-values obtained from Wilcoxon rank-sum tests
pwm_original = [pwm1, pwm2, pwm3, pwm4, pwm5, pwm6];

adjusted_p_values_WRS_Mothers = benjamini_hochberg_correction(pwm_original);

% % Identify significant comparisons
significant_comparisons_m = find(adjusted_p_values_WRS_Mothers <= FDR_threshold);


% Display the original p-values and their corrected values
fprintf('WRS-test Mothers \n');
fprintf('Comparison\tOriginal p-values\tBH-Corrected p-values\n');
for i = 1:length(pwm_original)
    fprintf('%s\t%f\t\t\t%f\n',group_comparisons{i}, pwm_original(i), adjusted_p_values_WRS_Mothers(i));
end

 
% Display significant comparisons
if ~isempty(significant_comparisons_m)
    fprintf('\nSignificant Comparisons (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_comparisons_m)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_comparisons_m(i), adjusted_p_values_WRS_Mothers(significant_comparisons_m(i)));
    end
else
    disp('No significant differences were found after FDR correction.');
end

%% boxplot infants only

XXc=[MM1c.HI_Ct;MM2c.HI_Ct;MM3c.HI_Ct;MM5c.HI_Ct];
YYc=[ones(size(MM1c.HI_Ct,1),1);2*ones(size(MM1c.HI_Ct,1),1);3*ones(size(MM1c.HI_Ct,1),1);4*ones(size(MM1c.HI_Ct,1),1)];


cols=1;
rows=1;
figure
xSize = 15; Xs=xSize; ySize = 13.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

boxplot(XXc,YYc,'Color','k','Symbol','');
set(gca,'Fontsize',60);box on;
hold on;
scatter(1*ones(1,length(MM1c.HI_Ct)),MM1c.HI_Ct,dot_size,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(1,length(MM2c.HI_Ct)),MM2c.HI_Ct,dot_size,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(3*ones(1,length(MM3c.HI_Ct)),MM3c.HI_Ct,dot_size,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(4*ones(1,length(MM5c.HI_Ct)),MM5c.HI_Ct,dot_size,'filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot([1,2,3,4],meanMHI.mean_HI_Ct,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20);
set(gca,'Fontsize',40);box on;
xticks([1,2,3,4])
ylim([18,45])
set(gca, 'XTickLabel', {'HI$\rightarrow$RSV','HI\&RSV','RSV$\rightarrow$HI','HI'},'TickLabelInterpreter','latex','Fontsize',25)
ylabel('HI Ct','Interpreter','latex')
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(30)

ylim([18,50])
set(gca,'Fontsize',60);box on;


if yxis_reverse==1
    set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
end

%%  t-test among infant groups
[h1c,p1c]=ttest2(MM1c.HI_Ct,MM2c.HI_Ct);
[h2c,p2c]=ttest2(MM1c.HI_Ct,MM3c.HI_Ct);
[h3c,p3c]=ttest2(MM1c.HI_Ct,MM5c.HI_Ct);
[h4c,p4c]=ttest2(MM2c.HI_Ct,MM3c.HI_Ct);
[h5c,p5c]=ttest2(MM2c.HI_Ct,MM5c.HI_Ct);
[h6c,p6c]=ttest2(MM3c.HI_Ct,MM5c.HI_Ct);


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
[pwc1,hwc1,statswc1] = ranksum(MM1c.HI_Ct,MM2c.HI_Ct);
[pwc2,hwc2,statswc2] = ranksum(MM1c.HI_Ct,MM3c.HI_Ct); 
[pwc3,hwc3,statswc3] = ranksum(MM1c.HI_Ct,MM5c.HI_Ct);
[pwc4,hwc4,statswc4] = ranksum(MM2c.HI_Ct,MM3c.HI_Ct);
[pwc5,hwc5,statswc5] = ranksum(MM2c.HI_Ct,MM5c.HI_Ct);
[pwc6,hwc6,statswc6] = ranksum(MM3c.HI_Ct,MM5c.HI_Ct); 

FDR_threshold=0.05;

% % Benjamini-Hochberg (False Discovery Rate, FDR) correction 

% Store the original p-values and results in a cell array
pwc_original = [pwc1, pwc2, pwc3, pwc4, pwc5, pwc6];

adjusted_p_values_WRS_Infants = benjamini_hochberg_correction(pwc_original);

% % Identify significant comparisons
significant_comparisons_c = find(adjusted_p_values_WRS_Infants <= FDR_threshold);


% Display the original p-values and their corrected values
fprintf('WRS-test Infants\n');
fprintf('Comparison\tOriginal p-values\tCorrected p-values\n');
for i = 1:length(pwc_original)
    fprintf('%s\t%f\t\t\t%f\n',group_comparisons{i}, pwc_original(i), adjusted_p_values_WRS_Infants(i));
end
 
% Display significant comparisons
if ~isempty(significant_comparisons_c)
    fprintf('\nSignificant Comparisons (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_comparisons_c)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_comparisons_c(i), adjusted_p_values_WRS_Infants(significant_comparisons_c(i)));
    end
else
    disp('No significant differences were found after FDR correction.');
end




