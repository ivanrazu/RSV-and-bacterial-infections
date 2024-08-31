clear
clc

dot_size=400;


% loads info from RSV->SP
load('MtxGroup3_child_with_demographic_data.mat','MtxGroup3_child');

M3=MtxGroup3_child;

% % Subjects who came up SP+ in the following visit immedially after detection of RSV
M3subjects1=[14;19;29;84;146;330;349;671;1810];

% %Subjects who came up SP+ in the second  or later visit after detection of RSV
M3subjects2=[31;225;259;344;352;411;1656];


for j=1:length(M3subjects1)
    indx=find(M3.subject_id==M3subjects1(j));
    A=M3.SP_Ct_Mean(indx);
    B=A(A>0);
    
    M3_SP_day1a(j,:)=B(1);
    if size(B,1)>1
        M3_SP_day2a(j,:)=B(2);
        if size(B,1)>2
            M3_SP_day3a(j,:)=B(3);
            
            if size(B,1)>3
                M3_SP_day4a(j,:)=B(4);
                
                if size(B,1)>4
                    M3_SP_day5a(j,:)=B(5);
                    
                    if size(B,1)>5
                        M3_SP_day6a(j,:)=B(6);
                    end
                    
                end
            end
            
        end
    end
    
end



M3_SP_day1a=nonzeros(M3_SP_day1a);
M3_SP_day2a=nonzeros(M3_SP_day2a);
M3_SP_day3a=nonzeros(M3_SP_day3a);
M3_SP_day4a=nonzeros(M3_SP_day4a);
M3_SP_day5a=nonzeros(M3_SP_day5a);
% M3_SP_day6a=nonzeros(M3_SP_day6a);

Ma=[M3_SP_day1a;M3_SP_day2a;M3_SP_day3a;M3_SP_day4a;M3_SP_day5a];



for j=1:length(M3subjects2)
    indx=find(M3.subject_id==M3subjects2(j));
    A=M3.SP_Ct_Mean(indx);
    B=A(A>0);
    
    M3_SP_day1b(j,:)=B(1);
    if size(B,1)>1
        M3_SP_day2b(j,:)=B(2);
        if size(B,1)>2
            M3_SP_day3b(j,:)=B(3);
            
            if size(B,1)>3
                M3_SP_day4b(j,:)=B(4);
                
                if size(B,1)>4
                    M3_SP_day5b(j,:)=B(5);
                    
                    if size(B,1)>5
                        M3_SP_day6b(j,:)=B(6);
                    end
                    
                end
            end
            
        end
    end
    
end



M3_SP_day1b=nonzeros(M3_SP_day1b);
M3_SP_day2b=nonzeros(M3_SP_day2b);
M3_SP_day3b=nonzeros(M3_SP_day3b);
M3_SP_day4b=nonzeros(M3_SP_day4b);
% M3_SP_day5b=nonzeros(M3_SP_day5b);

MSP_alldays=M3.SP_Ct_Mean;
MSP_alldays=MSP_alldays(MSP_alldays>0);


Mb=[M3_SP_day1b;M3_SP_day2b;M3_SP_day3b;M3_SP_day4b];



%%
figure
xSize = 15; Xs=xSize; ySize = 11.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[350 150 xSize*50 ySize*50]);

boxplot([  Ma; Mb],[ones(length(Ma),1);...
    2*ones(length(Mb),1)],'Color','k','Symbol','')
hold on
scatter(1*ones(1,length(Ma)),Ma,dot_size,'filled','b','MarkerFaceColor','k','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
scatter(2*ones(1,length(Mb)),Mb,dot_size,'filled','b','MarkerFaceColor','k','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
hold on
plot(1,mean(Ma),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
plot(2,mean(Mb),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',20)
hold on
set(gca,'Fontsize',60);box on;

ylabel('SP Ct','interpreter','latex')
set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)

yxis_reverse=1;

% % plots line to denote significant difference
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.1)])
xt=[1,2];

if yxis_reverse==0
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.14, 'k-','LineWidth',6)
    dim = [.52 .57 .9 .3];
    ylim([18,45])
end

if yxis_reverse==1
    set(gca, 'YDir', 'reverse')
    ylim([15,45])
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*0.56, 'k-','LineWidth',6)
    hold on
    yticks([20,30,40])
    set(gca, 'YTickLabel', {'20','30','40'})
    dim = [.52 .60 .9 .3];
end

name={'$ * $'};
annotation('textbox',dim,'String',name,'interpreter','latex','Fontsize',55,'Color', 'k','EdgeColor','none');

ax.YAxis.FontSize = 60;
ax.YLabel.FontSize = 60;


row1 = {'SP+ in' 'SP+ after'};
row2 = {'next visit'  '2 or more visits'};
row3 = {'post-RSV'  'post-RSV' };
labelArray = [row1; row2; row3]; 
labelArray = strjust(pad(labelArray),'left');

tickLabels = strtrim(sprintf('%s\\newline%s\\newline%s\n', labelArray{:}));
ax = gca(); 
ax.XTick = 1:2; 
% ax.XLim = [0,3];
ax.XTickLabel = tickLabels; 
ax.TickLabelInterpreter = 'tex';

ax.XAxis.FontSize =15;

xtickangle(0)

if yxis_reverse==0
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String',...
    'RSV$\rightarrow$SP','interpreter','latex','EdgeColor','none','FontSize',50)
end

if yxis_reverse==1
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String',...
    'RSV$\rightarrow$SP','interpreter','latex','EdgeColor','none','FontSize',50)
end


%saveas(gca,'boxplot_infants_SP_rightafter_RSV_v4.png')


%% Statistics

% t-test 
[h,p]=ttest2(Ma,Mb); %%% p=0.01 

% Wilcoxon rank sum test
[p1a,h1a]=ranksum(Ma,Mb); %%% p=0.0085


