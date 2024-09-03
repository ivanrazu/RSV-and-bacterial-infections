
% % Not enough data to make this analysis here

clear
clc

dot_size=400;


load('MtxGroup3_mother_with_demographic_data.mat','MtxGroup3_mother');

%% 
M3=MtxGroup3_mother;

M3subjects1=[];
M3subjects2=[49;219;342;343;344];

if ~isempty(M3subjects1)

for j=1:length(M3subjects1)
    indx=find(M3.subject_id==M3subjects1(j));
A=M3.MC_Ct_Mean(indx);
B=A(A>0);

M3_MC_day1a(j,:)=B(1);
if size(B,1)>1
    M3_MC_day2a(j,:)=B(2);
    if size(B,1)>2
        M3_MC_day3a(j,:)=B(3);

        if size(B,1)>3
            M3_MC_day4a(j,:)=B(4);

            if size(B,1)>4
                M3_MC_day5a(j,:)=B(5);

                if size(B,1)>5
                    M3_MC_day6a(j,:)=B(6);
                end

            end
        end

    end
end

end

M3_MC_day1a=nonzeros(M3_MC_day1a);
M3_MC_day2a=nonzeros(M3_MC_day2a);


else

M3_MC_day1a=[];
M3_MC_day2a=[];

end


Ma=[M3_MC_day1a;M3_MC_day2a];

% Work arround to be able to plot boxplot when MM1a_new is empty
if isempty(M3subjects1)
    Ma=[zeros(3,1)];   
end




%%

for j=1:length(M3subjects2)
    indx=find(M3.subject_id==M3subjects2(j));
A=M3.MC_Ct_Mean(indx);
B=A(A>0);

M3_MC_day1b(j,:)=B(1);
if size(B,1)>1
    M3_MC_day2b(j,:)=B(2);
    if size(B,1)>2
        M3_MC_day3b(j,:)=B(3);

        if size(B,1)>3
            M3_MC_day4b(j,:)=B(4);

            if size(B,1)>4
                M3_MC_day5b(j,:)=B(5);

                if size(B,1)>5
                    M3_MC_day6b(j,:)=B(6);
                end

            end
        end

    end
end

end

M3_MC_day1b=nonzeros(M3_MC_day1b);
M3_MC_day2b=nonzeros(M3_MC_day2b);
MMC_alldays=M3.MC_Ct_Mean;
MMC_alldays=MMC_alldays(MMC_alldays>0);


Mb=[M3_MC_day1b;M3_MC_day2b];

%%
col_bl=[0 0.4470 0.7410];

figure
xSize = 15; Xs=xSize; ySize = 11.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

boxplot([ Ma; Mb],[ones(length(Ma),1);...
    2*ones(length(Mb),1)],'Color','k','Symbol','')
hold on
scatter(1*ones(1,length(Ma)),Ma,dot_size,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(1,length(Mb)),Mb,dot_size,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot(1,mean(Ma),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
plot(2,mean(Mb),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
set(gca,'Fontsize',60);box on;
ylabel('MC Ct','interpreter','latex')
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)

row1 = {'MC+ in' 'MC+ after'};
row2 = {'next visit'  '2 or more visits'};
row3 = {'post-RSV'  'post-RSV' };
labelArray = [row1; row2; row3]; 
labelArray = strjust(pad(labelArray),'left');

tickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\n', labelArray{:}));
ax = gca(); 
ax.XTick = 1:2; 
ax.XTickLabel = tickLabels; 
ax.TickLabelInterpreter = 'tex';
ax.XAxis.FontSize =15;

 yxis_reverse=1;
if yxis_reverse==0
    ylim([18,45])
end

xtickangle(0)
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String',...
    'RSV$\rightarrow$MC','interpreter','latex','EdgeColor','none','FontSize',50)

if yxis_reverse==1
    set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
end


%% Statistics

% t-test 
[h,p]=ttest2(Ma,Mb);

% Wilcoxon rank sum test
[pwc, hwc, statswc] = ranksum(Ma,Mb);


