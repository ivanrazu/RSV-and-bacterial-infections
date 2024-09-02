%Plots boxplots of RSV CT on Mothers only, infants only, and both combined
clear
clc

warning off


load('MtxGroup1_child_with_demographic_data.mat','MtxGroup1_child');
load('MtxGroup2_child_with_demographic_data.mat','MtxGroup2_child');
load('MtxGroup3_child_with_demographic_data.mat','MtxGroup3_child');
load('MtxGroup4_child_with_demographic_data.mat','MtxGroup4_child');
load('MtxGroup5_child_with_demographic_data.mat','MtxGroup5_child');


load('MtxGroup1_mother_with_demographic_data.mat','MtxGroup1_mother');
% load('MtxGroup1_mother_with_demographic_data.mat','MtxGroup2_mother');
load('MtxGroup3_mother_with_demographic_data.mat','MtxGroup3_mother');
load('MtxGroup4_mother_with_demographic_data.mat','MtxGroup4_mother');
load('MtxGroup5_mother_with_demographic_data.mat','MtxGroup5_mother');


M1m=MtxGroup1_mother;
% M2m=MtxGroup2_mother;
M3m=MtxGroup3_mother;
M4m=MtxGroup4_mother;
M5m=MtxGroup5_mother;


M1c=MtxGroup1_child;
M2c=MtxGroup2_child;
M3c=MtxGroup3_child;
M4c=MtxGroup4_child;
M5c=MtxGroup5_child;


%%

RSV_positive_only_M1c=M1c.RSV_CT(M1c.RSV_CT<99);
mean_RSV_positive_M1c=mean(RSV_positive_only_M1c);

RSV_positive_only_M2c=M2c.RSV_CT(M2c.RSV_CT<99);
mean_RSV_positive_M2c=mean(RSV_positive_only_M2c);

RSV_positive_only_M3c=M3c.RSV_CT(M3c.RSV_CT<99);
mean_RSV_positive_M3c=mean(RSV_positive_only_M3c);

RSV_positive_only_M4c=M4c.RSV_CT(M4c.RSV_CT<99);
mean_RSV_positive_M4c=mean(RSV_positive_only_M4c);


RSV_positive_only_M1m=M1m.RSV_CT(M1m.RSV_CT<99);
mean_RSV_positive_M1m=mean(RSV_positive_only_M1m);

% RSV_positive_only_M2m=M2m.RSV_CT(M2m.RSV_CT<99);
% mean_RSV_positive_M2m=mean(RSV_positive_only_M2m);

RSV_positive_only_M3m=M3m.RSV_CT(M3m.RSV_CT<99);
mean_RSV_positive_M3m=mean(RSV_positive_only_M3m);

RSV_positive_only_M4m=M4m.RSV_CT(M4m.RSV_CT<99);
mean_RSV_positive_M4m=mean(RSV_positive_only_M4m);
    

MM1c=table(RSV_positive_only_M1c, ones(1,length(RSV_positive_only_M1c))');
MM2c=table(RSV_positive_only_M2c,2*ones(1,length(RSV_positive_only_M2c))');
MM3c=table(RSV_positive_only_M3c,3*ones(1,length(RSV_positive_only_M3c))');
MM4c=table(RSV_positive_only_M4c,4*ones(1,length(RSV_positive_only_M4c))');

MM1c=renamevars(MM1c,["Var2","RSV_positive_only_M1c"],["group","RSV_Ct"]);
MM2c=renamevars(MM2c,["Var2","RSV_positive_only_M2c"],["group","RSV_Ct"]);
MM3c=renamevars(MM3c,["Var2","RSV_positive_only_M3c"],["group","RSV_Ct"]);
MM4c=renamevars(MM4c,["Var2","RSV_positive_only_M4c"],["group","RSV_Ct"]);


MM1m=table(RSV_positive_only_M1m, ones(1,length(RSV_positive_only_M1m))');
% MM2m=table(RSV_positive_only_M2m,2*ones(1,length(RSV_positive_only_M2m))');
MM3m=table(RSV_positive_only_M3m,3*ones(1,length(RSV_positive_only_M3m))');
MM4m=table(RSV_positive_only_M4m,4*ones(1,length(RSV_positive_only_M4m))');

MM1m=renamevars(MM1m,["Var2","RSV_positive_only_M1m"],["group","RSV_Ct"]);
% MM2m=renamevars(MM2m,["Var2","RSV_positive_only_M2m"],["group","RSV_Ct"]);
MM3m=renamevars(MM3m,["Var2","RSV_positive_only_M3m"],["group","RSV_Ct"]);
MM4m=renamevars(MM4m,["Var2","RSV_positive_only_M4m"],["group","RSV_Ct"]);


MMm=[MM1m;MM3m;MM4m];
meanMMm = groupsummary(MMm,'group','mean');

MMc=[MM1c;MM3c;MM4c];
meanMMc = groupsummary(MMc,'group','mean');


col_bl=[0 0.4470 0.7410];
col_blk=[0,0,0];

 alpha = 0.05;  % Set significance level


%% boxplot mothers only

figure
xSize = 15; Xs=xSize; ySize = 11.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

h=boxplot([MM1m.RSV_Ct;MM3m.RSV_Ct;MM4m.RSV_Ct],...
    [ones(size(MM1m.RSV_Ct,1),1);...
    2*ones(size(MM3m.RSV_Ct,1),1);3*ones(size(MM4m.RSV_Ct,1),1)],'Color','k','Symbol','');

hold on;
scatter(1,MM1m.RSV_Ct,1000,'d','filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
% scatter(2,MM2m.RSV_Ct,1000,'filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% hold on
scatter(2,MM3m.RSV_Ct,1000,'d','filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on

scatter(3,MM4m.RSV_Ct,1000,'d','filled','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot([1,2,3],[mean(MM1m.RSV_Ct(MM1m.RSV_Ct>0)),...
    mean(MM3m.RSV_Ct(MM3m.RSV_Ct>0)),mean(MM4m.RSV_Ct(MM4m.RSV_Ct>0))],...
    'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20);
hold on

yxis_reverse=1;

if yxis_reverse==0
    ylim([16,45])
end

set(gca,'Fontsize',35);box on;

set(gca, 'XTickLabel', {'SA$\rightarrow$RSV','RSV$\rightarrow$SA','RSV'},'TickLabelInterpreter','latex','Fontsize',50)
ylabel('RSV Ct','Interpreter','latex')
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)

xtickangle(30)

if yxis_reverse==1
set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
end


set(gca,'Fontsize',60);box on;


 %%   T-test among mother groups

[h2m,p2m]=ttest2(MM1m.RSV_Ct,MM3m.RSV_Ct);
[h3m,p3m]=ttest2(MM1m.RSV_Ct,MM4m.RSV_Ct);
[h6m,p6m]=ttest2(MM3m.RSV_Ct,MM4m.RSV_Ct);

% p-values obtained from t-test
pim=[p2m,p3m,p6m];

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

group_comparisons={ 'SA->RSV	vs RSV->SA ' , 'SA->RSV vs RSV   ',   'RSV->SA vs RSV   '};


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

%% Wilcoxon rank sum test in Mothers
[pwm2,hwm2,statswm2] = ranksum(MM1m.RSV_Ct,MM3m.RSV_Ct);
[pwm3,hwm3,statswm3] = ranksum(MM1m.RSV_Ct,MM4m.RSV_Ct); 
[pwm6,hwm6,statswm6] = ranksum(MM3m.RSV_Ct,MM4m.RSV_Ct); 

FDR_threshold=0.05;

% Benjamini-Hochberg (False Discovery Rate, FDR) correction 

% p-values obtained from Wilcoxon rank-sum tests
pwm_original = [pwm2, pwm3, pwm6];

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
    fprintf('\nSignificant Comparisons according to Wilcoxon rank sum test (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_indices_WRS_test_BH_m)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_indices_WRS_test_BH_m(i), adjusted_p_values_WRS_Mothers(significant_indices_WRS_test_BH_m(i)));
    end
else
    disp('No significant differences were found according to Wilcoxon rank sum test with FDR correction.');
end

%% boxplot infants only

MMc=[MM1c;MM2c;MM3c;MM4c];
meanMMc = groupsummary(MMc,'group','mean');

figure
xSize = 15; Xs=xSize; ySize = 11.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

h=boxplot([MM1c.RSV_Ct;MM2c.RSV_Ct;MM3c.RSV_Ct;MM4c.RSV_Ct],...
    [ones(size(MM1c.RSV_Ct,1),1);2*ones(size(MM2c.RSV_Ct,1),1);...
    3*ones(size(MM3c.RSV_Ct,1),1);4*ones(size(MM4c.RSV_Ct,1),1)],'Color','k','Symbol','');

hold on;
scatter(1,MM1c.RSV_Ct,1000,'d','filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2,MM2c.RSV_Ct,1000,'d','filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(3,MM3c.RSV_Ct,1000, 'd','filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(4,MM4c.RSV_Ct,1000, 'd','filled','MarkerFaceColor',col_blk,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot([1,2,3,4],[mean(MM1c.RSV_Ct(MM1c.RSV_Ct>0)),mean(MM2c.RSV_Ct(MM2c.RSV_Ct>0)),...
    mean(MM3c.RSV_Ct(MM3c.RSV_Ct>0)),mean(MM4c.RSV_Ct(MM4c.RSV_Ct>0))],...
    'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20);

hold on
if yxis_reverse==0
    ylim([16,45])
end

set(gca,'Fontsize',35);box on;

set(gca, 'XTickLabel', {'SA$\rightarrow$RSV','SA\&RSV','RSV$\rightarrow$SA','RSV'},'TickLabelInterpreter','latex','FontSize',50)

ylabel('RSV Ct','Interpreter','latex')
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)
xtickangle(30)

if yxis_reverse==1
set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
end

set(gca,'Fontsize',60);box on;

%% t-test among infant groups
[h1c,p1c]=ttest2(MM1c.RSV_Ct,MM3c.RSV_Ct);
[h2c,p2c]=ttest2(MM1c.RSV_Ct,MM4c.RSV_Ct);
[h3c,p3c]=ttest2(MM2c.RSV_Ct,MM3c.RSV_Ct);
[h4c,p4c]=ttest2(MM2c.RSV_Ct,MM4c.RSV_Ct);
[h5c,p5c]=ttest2(MM3c.RSV_Ct,MM4c.RSV_Ct);

pic=[p1c,p2c,p3c,p4c,p5c];

% Find which p-values are lower than alpha
significant_indices_ttest_c = find(pic < alpha);

% Create a cell array of variable names
variable_names_c = {'p1c', 'p2c', 'p3c', 'p4c', 'p5c'};

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
pct_original = [p1c, p2c, p3c, p4c, p5c];

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

%% %% Wilcoxon rank sum test among infants 
[pwc1,hwc1,statswc1] = ranksum(MM1c.RSV_Ct,MM2c.RSV_Ct);
[pwc2,hwc2,statswc2] = ranksum(MM1c.RSV_Ct,MM3c.RSV_Ct);
[pwc3,hwc3,statswc3] = ranksum(MM1c.RSV_Ct,MM4c.RSV_Ct);
[pwc4,hwc4,statswc4] = ranksum(MM2c.RSV_Ct,MM3c.RSV_Ct);
[pwc5,hwc5,statswc5] = ranksum(MM2c.RSV_Ct,MM4c.RSV_Ct);
[pwc6,hwc6,statswc6] = ranksum(MM3c.RSV_Ct,MM4c.RSV_Ct); 

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
    fprintf('\nSignificant Comparisons according to Wilcoxon rank sum test (FDR-corrected) at a threshold of %f:\n', FDR_threshold);
    for i = 1:length(significant_indices_WRS_BH_test_c)
        fprintf('Comparison %d (p = %f) is significant.\n', significant_indices_WRS_BH_test_c(i), adjusted_p_values_WRS_Infants(significant_indices_WRS_BH_test_c(i)));
    end
else
    disp('No significant differences were found according to Wilcoxon rank sum test with FDR correction.');
end


