
%% Not enough data to make this analysis here

% % clear
% % clc
% % 
% % load('MtxGroup3_mother.mat','MtxGroup3_mother');
% % 
% % 
% % M3=MtxGroup3_mother;
% % 
% % M3subjects1=[39;219;1665];
% % M3subjects2=[49;71;134;180;342;719;1817];
% % 
% % 
% % 
% % for j=1:length(M3subjects1)
% %     indx=find(M3.subject_id==M3subjects1(j));
% % A=M3.SP_Ct_Mean(indx);
% % B=A(A>0);
% % 
% % M3_SP_day1a(j,:)=B(1);
% % if size(B,1)>1
% %     M3_SP_day2a(j,:)=B(2);
% %     if size(B,1)>2
% %         M3_SP_day3a(j,:)=B(3);
% %         
% %         if size(B,1)>3
% %             M3_SP_day4a(j,:)=B(4);
% %             
% %             if size(B,1)>4
% %                 M3_SP_day5a(j,:)=B(5);
% %                 
% %                 if size(B,1)>5
% %                     M3_SP_day6a(j,:)=B(6);
% %                 end
% %                 
% %             end
% %         end
% %         
% %     end
% % end
% % 
% % end
% % 
% % 
% % 
% % M3_SP_day1a=nonzeros(M3_SP_day1a);
% % M3_SP_day2a=nonzeros(M3_SP_day2a);
% % 
% % Ma=[M3_SP_day1a;M3_SP_day2a];
% % 
% % for j=1:length(M3subjects2)
% %     indx=find(M3.subject_id==M3subjects2(j));
% % A=M3.SP_Ct_Mean(indx);
% % B=A(A>0);
% % 
% % M3_SP_day1b(j,:)=B(1);
% % if size(B,1)>1
% %     M3_SP_day2b(j,:)=B(2);
% %     if size(B,1)>2
% %         M3_SP_day3b(j,:)=B(3);
% %         
% %         if size(B,1)>3
% %             M3_SP_day4b(j,:)=B(4);
% %             
% %             if size(B,1)>4
% %                 M3_SP_day5b(j,:)=B(5);
% %                 
% %                 if size(B,1)>5
% %                     M3_SP_day6b(j,:)=B(6);
% %                 end
% %                 
% %             end
% %         end
% %         
% %     end
% % end
% % 
% % end
% % 
% % 
% % 
% % M3_SP_day1b=nonzeros(M3_SP_day1b);
% % M3_SP_day2b=nonzeros(M3_SP_day2b);
% % M3_SP_day3b=nonzeros(M3_SP_day3b);
% % 
% % MSP_alldays=M3.SP_Ct_Mean;
% % MSP_alldays=MSP_alldays(MSP_alldays>0);
% % 
% % 
% % Mb=[M3_SP_day1b;M3_SP_day2b;M3_SP_day3b];
% % 
% % %%
% % col_bl=[0 0.4470 0.7410];
% % 
% % figure
% % xSize = 15; Xs=xSize; ySize = 11.5;xLeft = (xSize-xSize)/2; Ys=ySize; yTop = (ySize-ySize)/2;
% % set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
% % set(gcf,'Position',[350 150 xSize*50 ySize*50]);
% % 
% % boxplot([ Ma; Mb],[ones(length(Ma),1);...
% %     2*ones(length(Mb),1)],'Color','k','Symbol','')
% % hold on
% % scatter(1*ones(1,length(Ma)),Ma,200,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% % hold on
% % scatter(2*ones(1,length(Mb)),Mb,200,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% % hold on
% % plot(1,mean(Ma),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
% % hold on
% % plot(2,mean(Mb),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
% % ylim([18,40])
% % set(gca,'Fontsize',40);box on;
% % 
% % % xticks((1:1:5))
% % % set(gca, 'XTickLabel', {'RSV->SP visit 1','RSV->SP visit 2 and later'})
% % % set(gca, 'XTickLabel', {'a','b','c','d','e'})
% % % set(gca, 'XTickLabel', {'RSV$\rightarrow$SP','(a)','(b)','(c)','(d)'},'TickLabelInterpreter','latex')
% % 
% % 
% % ylabel('SP Ct','interpreter','latex')
% % set(findobj(gca,'type','line'),'linew',4)
% % set(gca,'linew',4)
% % 
% % 
% % return
% % 
% % %%
% % figure
% % 
% % boxplot([ [Ma;Mb]; Ma; Mb ;M3_SP_day1a; M3_SP_day1b],[ones(length([Ma;Mb]),1);...
% %     2*ones(length(Ma),1);3*ones(length(Mb),1);4*ones(length(M3_SP_day1a),1);5*ones(length(M3_SP_day1b),1)],'Color','k','Symbol','')
% % hold on
% % scatter(ones(1,length(Ma)),Ma,200,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% % hold on
% % scatter(ones(1,length(Mb)),Mb,200,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% % hold on
% % scatter(2*ones(1,length(Ma)),Ma,200,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% % hold on
% % scatter(3*ones(1,length(Mb)),Mb,200,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% % hold on
% % scatter(4*ones(1,length(M3_SP_day1a)),M3_SP_day1a,200,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% % hold on
% % scatter(5*ones(1,length(M3_SP_day1b)),M3_SP_day1b,200,'filled','b','MarkerFaceColor',col_bl,'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);
% % hold on
% % plot(1,mean([Ma;Mb]),'r*','LineWidth',6)
% % hold on
% % plot(2,mean(Ma),'r*','LineWidth',6)
% % hold on
% % plot(3,mean(Mb),'r*','LineWidth',6)
% % hold on
% % plot(4,mean(M3_SP_day1a),'r*','LineWidth',6)
% % hold on
% % plot(5,mean(M3_SP_day1b),'r*','LineWidth',6)
% % hold on
% % xline(1.5,'k-','LineWidth',1')
% % hold on
% % xline(3.5,'k-','LineWidth',1')
% % ylim([18,40])
% % set(gca,'Fontsize',40);box on;
% % xticks((1:1:5))
% % % set(gca, 'XTickLabel', {'RSV->SP visit 1','RSV->SP visit 2 and later'})
% % % set(gca, 'XTickLabel', {'a','b','c','d','e'})
% % % set(gca, 'XTickLabel', {'RSV$\rightarrow$SP','(a)','(b)','(c)','(d)'},'TickLabelInterpreter','latex')
% % 
% % ylabel('SP Ct','interpreter','latex')
% % ylabel('SP Ct','interpreter','latex')
% % set(findobj(gca,'type','line'),'linew',4)
% % set(gca,'linew',4)
% % [h1,p1]=ttest2(Ma,Mb)
% % [h,p]=ttest2(M3_SP_day1a,M3_SP_day1b)
% % 
