clear 
clc
load('groupid_info_child_with_demograp_data.mat','subjects_child','group1id_child','group2id_child',...
    'group3id_child','group4id_child','group5id_child','group6id_child','Mtx_new_child');

dot_size=200;


% Transforms cell type into mat  
Mtx_new_child.Abx_last_visit=cell2mat(Mtx_new_child.Abx_last_visit);
Mtx_new_child.Abx_curr_visit=cell2mat(Mtx_new_child.Abx_curr_visit);

group=group1id_child;
MtxGroup1_child=[];
for j=1:length(group)
[indx,pos]=find(Mtx_new_child.subject_id==group(j));
MtxGroup1_child=[Mtx_new_child(indx,:);MtxGroup1_child];
end

MtxGroup1_child=[MtxGroup1_child, table(ones(height(MtxGroup1_child),1))];
MtxGroup1_child= renamevars(MtxGroup1_child,'Var1','Group');

 
%%

group=group2id_child;
MtxGroup2_child=[];
for j=1:length(group)
[indx,pos]=find(Mtx_new_child.subject_id==group(j));
MtxGroup2_child=[Mtx_new_child(indx,:);MtxGroup2_child];

end

MtxGroup2_child=[MtxGroup2_child, table(2*ones(height(MtxGroup2_child),1))];
MtxGroup2_child= renamevars(MtxGroup2_child,'Var1','Group');

%%

group=group3id_child;
MtxGroup3_child=[];
for j=1:length(group)
[indx,pos]=find(Mtx_new_child.subject_id==group(j));
MtxGroup3_child=[Mtx_new_child(indx,:);MtxGroup3_child];

end

MtxGroup3_child=[MtxGroup3_child, table(3*ones(height(MtxGroup3_child),1))];
MtxGroup3_child= renamevars(MtxGroup3_child,'Var1','Group');


%%

group=group4id_child;
MtxGroup4_child=[];
for j=1:length(group)
[indx,pos]=find(Mtx_new_child.subject_id==group(j));
MtxGroup4_child=[Mtx_new_child(indx,:);MtxGroup4_child];

end

MtxGroup4_child=[MtxGroup4_child, table(4*ones(height(MtxGroup4_child),1))];
MtxGroup4_child= renamevars(MtxGroup4_child,'Var1','Group');


%%

group=group5id_child;
MtxGroup5_child=[];
for j=1:length(group)
[indx,pos]=find(Mtx_new_child.subject_id==group(j));
MtxGroup5_child=[Mtx_new_child(indx,:);MtxGroup5_child];

end

MtxGroup5_child=[MtxGroup5_child, table(5*ones(height(MtxGroup5_child),1))];
MtxGroup5_child= renamevars(MtxGroup5_child,'Var1','Group');


%%

group=group6id_child;
MtxGroup6_child=[];
for j=1:length(group)
[indx,pos]=find(Mtx_new_child.subject_id==group(j));
MtxGroup6_child=[Mtx_new_child(indx,:);MtxGroup6_child];

end

MtxGroup6_child=[MtxGroup6_child, table(6*ones(height(MtxGroup6_child),1))];
MtxGroup6_child= renamevars(MtxGroup6_child,'Var1','Group');


%%

% % All groups together Exlcluding group4 group6
Mxt_child_by_group = [MtxGroup1_child;MtxGroup2_child;MtxGroup3_child;MtxGroup4_child;MtxGroup5_child];

% Removes -1 in Abx data and replace them by 0
Mxt_child_by_group.Abx_last_visit(Mxt_child_by_group.Abx_last_visit==-1)=0;
Mxt_child_by_group.Abx_curr_visit(Mxt_child_by_group.Abx_curr_visit==-1)=0;

% Not including reported values at last visit as are inconsistent
Abx =  0*Mxt_child_by_group.Abx_last_visit + Mxt_child_by_group.Abx_curr_visit;
Mxt_child_by_group = [Mxt_child_by_group, table(Abx)];

M= [table(Mxt_child_by_group.subject_id,'VariableNames',{'subject_id'}), ...
    table(Mxt_child_by_group.Abx,'VariableNames',{'Abx'}),...
    table(Mxt_child_by_group.SP_Ct_Mean,'VariableNames',{'SP_CT'}),...
    table(Mxt_child_by_group.RSV_CT,'VariableNames',{'RSV_CT'}),...
    table(Mxt_child_by_group.Group,'VariableNames',{'Group'})];

idx_G1=find(M.Group==1);
idx_G2=find(M.Group==2);
idx_G3=find(M.Group==3);
idx_G4=find(M.Group==4);
idx_G5=find(M.Group==5);

M1= M(idx_G1,:);
M2= M(idx_G2,:);
M3= M(idx_G3,:);
M4= M(idx_G4,:);
M5= M(idx_G5,:);

%%   Matrix with SP CT of subjects that were given Abx 
M1a=[];
for j=1:length(group1id_child)    
    inx=find(M1.subject_id==group1id_child(j));    
    first_Abx=find(M1.Abx(inx,:)>0);    
    % Colecting CT values starting on the following visit after Abx administration
    SP_CT_first_Abx=M1.SP_CT(inx(first_Abx+1:end));    

    M1a = [ M1a;  [M1.subject_id(inx(first_Abx+1:end))  M1.Abx(inx(first_Abx+1:end))   SP_CT_first_Abx    ]] ;     
end

M1a_new = [table(M1a(:,1),'VariableNames',{'subject_id'}),...
           table(M1a(:,2),'VariableNames',{'Abx'}),...
           table(M1a(:,3),'VariableNames',{'SP_CT'})];
todelete_M1a_new = M1a_new.SP_CT==-1;
M1a_new(todelete_M1a_new,:)=[];
todelete2_M1a_new = M1a_new.SP_CT==0;
M1a_new(todelete2_M1a_new,:)=[];

%% Matrix with SP CT of subjects that were NOT given Abx 
s=1;
for  j=1:length(group1id_child)    
    in=find(M1.subject_id==group1id_child(j));    
    if sum(M1.Abx(in,:))==0
        sub_M1_no_Abx(s)=group1id_child(j);
        s=s+1;        
    end    
end

clear inx
clear first_Abx
clear SP_CT_first_Abx
M1b=[];
for j=1:length(sub_M1_no_Abx)    
    inx=find(M1.subject_id==sub_M1_no_Abx(j));    
    first_Abx=find(M1.Abx(inx,:)==0);    
    SP_CT_first_Abx=M1.SP_CT(inx(first_Abx));    
    M1b = [ M1b;  [M1.subject_id(inx(first_Abx))  M1.Abx(inx(first_Abx))   SP_CT_first_Abx  ]] ;
     
end

M1b_new = [table(M1b(:,1),'VariableNames',{'subject_id'}), table(M1b(:,2),'VariableNames',{'Abx'}), table(M1b(:,3),'VariableNames',{'SP_CT'})];
todelete_M1b_new = M1b_new.SP_CT==-1;
M1b_new(todelete_M1b_new,:)=[];
todelete2_M1b_new = M1b_new.SP_CT==0;
M1b_new(todelete2_M1b_new,:)=[];



%%  Matrix with SP CT of subjects that were given Abx 
M2a=[];
for j=1:length(group2id_child)    
    inx=find(M2.subject_id==group2id_child(j));    
    first_Abx=find(M2.Abx(inx,:)>0);    
    % Colecting CT values starting on the following visit after Abx administration
    SP_CT_first_Abx=M2.SP_CT(inx(first_Abx+1:end));    
    M2a = [ M2a;  [M2.subject_id(inx(first_Abx+1:end))  M2.Abx(inx(first_Abx+1:end))   SP_CT_first_Abx    ]] ;     
end

M2a_new = [table(M2a(:,1),'VariableNames',{'subject_id'}),...
           table(M2a(:,2),'VariableNames',{'Abx'}),...
           table(M2a(:,3),'VariableNames',{'SP_CT'})];
todelete_M2a_new = M2a_new.SP_CT==-1;
M2a_new(todelete_M2a_new,:)=[];
todelete2_M2a_new = M2a_new.SP_CT==0;
M2a_new(todelete2_M2a_new,:)=[];

%% Matrix with SP CT of subjects that were NOT given Abx 
s=1;
for  j=1:length(group2id_child)    
    in=find(M2.subject_id==group2id_child(j));    
    if sum(M2.Abx(in,:))==0
        sub_M2_no_Abx(s)=group2id_child(j);
        s=s+1;        
    end    
end

clear inx
clear first_Abx
clear SP_CT_first_Abx
M2b=[];
for j=1:length(sub_M2_no_Abx)    
    inx=find(M2.subject_id==sub_M2_no_Abx(j));    
    first_Abx=find(M2.Abx(inx,:)==0);    
    SP_CT_first_Abx=M2.SP_CT(inx(first_Abx));    
    M2b = [ M2b;  [M2.subject_id(inx(first_Abx))  M2.Abx(inx(first_Abx))   SP_CT_first_Abx  ]] ;
     
end

M2b_new = [table(M2b(:,1),'VariableNames',{'subject_id'}), table(M2b(:,2),'VariableNames',{'Abx'}), table(M2b(:,3),'VariableNames',{'SP_CT'})];
todelete_M2b_new = M2b_new.SP_CT==-1;
M2b_new(todelete_M2b_new,:)=[];
todelete2_M2b_new = M2b_new.SP_CT==0;
M2b_new(todelete2_M2b_new,:)=[];

%%  Matrix with SP CT of subjects that were given Abx 
M3a=[];

for j=1:length(group3id_child)    
    inx=find(M3.subject_id==group3id_child(j));    
    first_Abx=find(M3.Abx(inx,:)>0);    
    % Colecting CT values starting on the following visit after Abx administration   
    SP_CT_first_Abx=M3.SP_CT(inx(first_Abx+1:end));    
    M3a = [ M3a;  [M3.subject_id(inx(first_Abx+1:end))  M3.Abx(inx(first_Abx+1:end))   SP_CT_first_Abx    ]] ;     
end

M3a_new = [table(M3a(:,1),'VariableNames',{'subject_id'}),...
           table(M3a(:,2),'VariableNames',{'Abx'}),...
           table(M3a(:,3),'VariableNames',{'SP_CT'})];
todelete_M3a_new = M3a_new.SP_CT==-1;
M3a_new(todelete_M3a_new,:)=[];
todelete2_M3a_new = M3a_new.SP_CT==0;
M3a_new(todelete2_M3a_new,:)=[];

%% Matrix with SP CT of subjects that were NOT given Abx 
s=1;
for  j=1:length(group3id_child)    
    in=find(M3.subject_id==group3id_child(j));    
    if sum(M3.Abx(in,:))==0
        sub_M3_no_Abx(s)=group3id_child(j);
        s=s+1;        
    end    
end

clear inx
clear first_Abx
clear SP_CT_first_Abx
M3b=[];

for j=1:length(sub_M3_no_Abx)    
    inx=find(M3.subject_id==sub_M3_no_Abx(j));    
    first_Abx=find(M3.Abx(inx,:)==0);    
    SP_CT_first_Abx=M3.SP_CT(inx(first_Abx));    
    M3b = [ M3b;  [M3.subject_id(inx(first_Abx))  M3.Abx(inx(first_Abx))   SP_CT_first_Abx  ]] ;
     
end

M3b_new = [table(M3b(:,1),'VariableNames',{'subject_id'}), table(M3b(:,2),'VariableNames',{'Abx'}), table(M3b(:,3),'VariableNames',{'SP_CT'})];
todelete_M3b_new = M3b_new.SP_CT==-1;
M3b_new(todelete_M3b_new,:)=[];
todelete2_M3b_new = M3b_new.SP_CT==0;
M3b_new(todelete2_M3b_new,:)=[];


%% Matrix with SP CT of subjects that were given Abx 
M4a=[];
for j=1:length(group4id_child)    
    inx=find(M4.subject_id==group4id_child(j));    
    first_Abx=find(M4.Abx(inx,:)>0);    
    % Colecting CT values starting on the following visit after Abx administration
    SP_CT_first_Abx=M4.SP_CT(inx(first_Abx+1:end));    
    M4a = [ M4a;  [M4.subject_id(inx(first_Abx+1:end))  M4.Abx(inx(first_Abx+1:end))   SP_CT_first_Abx    ]] ;     
end

M4a_new = [table(M4a(:,1),'VariableNames',{'subject_id'}),...
           table(M4a(:,2),'VariableNames',{'Abx'}),...
           table(M4a(:,3),'VariableNames',{'SP_CT'})];
todelete_M4a_new = M4a_new.SP_CT==-1;
M4a_new(todelete_M4a_new,:)=[];
todelete2_M4a_new = M4a_new.SP_CT==0;
M4a_new(todelete2_M4a_new,:)=[];

%% Matrix with SP CT of subjects that were NOT given Abx 
s=1;
for  j=1:length(group4id_child)    
    in=find(M4.subject_id==group4id_child(j));    
    if sum(M4.Abx(in,:))==0
        sub_M4_no_Abx(s)=group4id_child(j);
        s=s+1;        
    end    
end

clear inx
clear first_Abx
clear SP_CT_first_Abx
M4b=[];
for j=1:length(sub_M4_no_Abx)    
    inx=find(M4.subject_id==sub_M4_no_Abx(j));    
    first_Abx=find(M4.Abx(inx,:)==0);    
    SP_CT_first_Abx=M4.SP_CT(inx(first_Abx));    
    M4b = [ M4b;  [M4.subject_id(inx(first_Abx))  M4.Abx(inx(first_Abx))   SP_CT_first_Abx  ]] ;
     
end

M4b_new = [table(M4b(:,1),'VariableNames',{'subject_id'}), table(M4b(:,2),'VariableNames',{'Abx'}), table(M4b(:,3),'VariableNames',{'SP_CT'})];
todelete_M4b_new = M4b_new.SP_CT==-1;
M4b_new(todelete_M4b_new,:)=[];
todelete2_M4b_new = M4b_new.SP_CT==0;
M4b_new(todelete2_M4b_new,:)=[];



%% Matrix with SP CT of subjects that were given Abx 
M5a=[];
clear inx
clear first_Abx
clear SP_CT_first_Abx
for j=1:length(group5id_child)    
    inx=find(M5.subject_id==group5id_child(j));    
    first_Abx=find(M5.Abx(inx,:)>0);  
    % Colecting CT values starting on the following visit after Abx administration    
    SP_CT_first_Abx=M5.SP_CT(inx(first_Abx+1:end));    
    M5a = [ M5a;  [M5.subject_id(inx(first_Abx+1:end))  M5.Abx(inx(first_Abx+1:end))   SP_CT_first_Abx  ]] ;      
end

M5a_new = [table(M5a(:,1),'VariableNames',{'subject_id'}), table(M5a(:,2),'VariableNames',{'Abx'}), table(M5a(:,3),'VariableNames',{'SP_CT'})];
todelete_M5a_new = M5a_new.SP_CT==-1;
M5a_new(todelete_M5a_new,:)=[];
todelete2_M5a_new = M5a_new.SP_CT==0;
M5a_new(todelete2_M5a_new,:)=[];

%% Matrix with SP CT of subjects that were NOT given Abx 

clear in
r=1;
for  j=1:length(group5id_child)    
    in=find(M5.subject_id==group5id_child(j));    
    if sum(M5.Abx(in,:))==0
        sub_M5_no_Abx(r)=group5id_child(j);
        r=r+1;        
    end    
end

clear inx
clear first_Abx
clear SP_CT_first_Abx
M5b=[];
for j=1:length(sub_M5_no_Abx)    
    inx=find(M5.subject_id==sub_M5_no_Abx(j));    
    first_Abx=find(M5.Abx(inx,:)==0);    
    SP_CT_first_Abx=M5.SP_CT(inx(first_Abx));    
    M5b = [ M5b;  [M5.subject_id(inx(first_Abx))  M5.Abx(inx(first_Abx))   SP_CT_first_Abx  ]] ;     
end

M5b_new = [table(M5b(:,1),'VariableNames',{'subject_id'}), table(M5b(:,2),'VariableNames',{'Abx'}), table(M5b(:,3),'VariableNames',{'SP_CT'})];
todelete_M5b_new = M5b_new.SP_CT==-1;
M5b_new(todelete_M5b_new,:)=[];
todelete2_M5b_new = M5b_new.SP_CT==0;
M5b_new(todelete2_M5b_new,:)=[];

%%

MM1 = [M1a_new.SP_CT; M1b_new.SP_CT];
MM2 = [M2a_new.SP_CT; M2b_new.SP_CT];
MM3 = [M3a_new.SP_CT; M3b_new.SP_CT];
MM4 = [M4a_new.SP_CT; M4b_new.SP_CT];
MM5= [ M5a_new.SP_CT; M5b_new.SP_CT];

VV1=[ones(length(M1a_new.SP_CT),1); 2*ones(length(M1b_new.SP_CT),1) ];
VV2=[ones(length(M2a_new.SP_CT),1); 2*ones(length(M2b_new.SP_CT),1) ];
VV3=[ones(length(M3a_new.SP_CT),1); 2*ones(length(M3b_new.SP_CT),1) ];
VV4=[ones(length(M4a_new.SP_CT),1); 2*ones(length(M4b_new.SP_CT),1) ];
VV5=[ones(length(M5a_new.SP_CT),1); 2*ones(length(M5b_new.SP_CT),1) ];

%% T-test

[h1,p1]=ttest2(M1a_new.SP_CT,M1b_new.SP_CT);%%% p= 0.1227
[h2,p2]=ttest2(M2a_new.SP_CT,M2b_new.SP_CT);%%% p=0.6578
[h3,p3]=ttest2(M3a_new.SP_CT,M3b_new.SP_CT);%%% p=0.7979
[h4,p4]=ttest2(M5a_new.SP_CT,M5b_new.SP_CT);%%% p=0.0259
%%

[p1a,h1a]=ranksum(M1a_new.SP_CT,M1b_new.SP_CT); %% p=0.113945
[p2a,h2a]=ranksum(M2a_new.SP_CT,M2b_new.SP_CT); %% p=0.785028
[p3a,h3a]=ranksum(M3a_new.SP_CT,M3b_new.SP_CT); %% p=0.770335
[p4a,h4a]=ranksum(M5a_new.SP_CT,M5b_new.SP_CT); %% p=0.024974

pmc=[p1a,p2a,p3a,p4a];

% Check for significant differences
alpha = 0.05;  % Set significance level
for i=1:4
    if pmc(i) < alpha
        fprintf('Infants Abx vs No Abx Group %d are significantly different p=%f.\n', i, pmc(i) );
    else
        fprintf('Infants Abx vs No Abx Group %d are NOT significantly different p=%f.\n', i, pmc(i) );
    end
end



%% T-test

[h1,p1]=ttest2(M3a_new.SP_CT,M3b_new.SP_CT)
[h2,p2]=ttest2(M5a_new.SP_CT,M5b_new.SP_CT)
[h3,p3]=ttest2(M3a_new.SP_CT,M5a_new.SP_CT)
[h4,p4]=ttest2(M3b_new.SP_CT,M5b_new.SP_CT)

%%
rng(123)

col_bl = [0 0.4470 0.7410];
col_blk = [0, 0, 0];

n1 = max([length(M1a_new.SP_CT), length(M2a_new.SP_CT), length(M3a_new.SP_CT), length(M5a_new.SP_CT)]);
n2 = max([length(M1b_new.SP_CT), length(M2b_new.SP_CT), length(M3b_new.SP_CT), length(M5b_new.SP_CT)]);

M1a_new.SP_CT(end+1:n1) = nan;
M2a_new.SP_CT(end+1:n1) = nan;
M3a_new.SP_CT(end+1:n1) = nan;
M5a_new.SP_CT(end+1:n1) = nan;

M1b_new.SP_CT(end+1:n2) = nan;
M2b_new.SP_CT(end+1:n2) = nan;
M3b_new.SP_CT(end+1:n2) = nan;
M5b_new.SP_CT(end+1:n2) = nan;

X = {[M1a_new.SP_CT, M2a_new.SP_CT, M3a_new.SP_CT, M5a_new.SP_CT], [M1b_new.SP_CT, M2b_new.SP_CT, M3b_new.SP_CT, M5b_new.SP_CT]};

meanM1235a = [mean(nonzeros(M1a_new.SP_CT), "omitnan"), mean(nonzeros(M2a_new.SP_CT), "omitnan"), mean(nonzeros(M3a_new.SP_CT), "omitnan"), mean(nonzeros(M5a_new.SP_CT), "omitnan")];
meanM1235b = [mean(nonzeros(M1b_new.SP_CT), "omitnan"), mean(nonzeros(M2b_new.SP_CT), "omitnan"), mean(nonzeros(M3b_new.SP_CT), "omitnan"), mean(nonzeros(M5b_new.SP_CT), "omitnan")];

cols = 1;
rows = 1;
fig = figure;
xSize = cols * 18;
Xs = xSize;
ySize = rows * 12;
xLeft = (xSize - xSize) / 2;
Ys = ySize;
yTop = (ySize - ySize) / 2;
set(gcf, 'PaperPosition', [xLeft yTop xSize ySize]);
set(gcf, 'Position', [100 25 xSize * 50 ySize * 55]);
boxplotGroup(X, 'Colors', [col_blk; col_blk; col_blk; col_blk], 'GroupType', 'betweenGroups', 'Symbol', '', 'groupLines', true);

jitterAmount = 0.15;
% dot_size = 300;

hold on;
scatter(1*ones(1,length(M1a_new.SP_CT)), M1a_new.SP_CT, dot_size, 'filled', 'MarkerFaceColor', col_blk, 'MarkerFaceAlpha', 0.6,'jitter','on','jitterAmount',0.15);
hold on;
scatter(2*ones(1,length(M1b_new.SP_CT)), M1b_new.SP_CT, dot_size, 'filled', 'MarkerFaceColor', col_blk, 'MarkerFaceAlpha', 0.6,'jitter','on','jitterAmount',0.15);
hold on;
scatter(4*ones(1,length(M2a_new.SP_CT)), M2a_new.SP_CT, dot_size, 'filled', 'MarkerFaceColor', col_blk, 'MarkerFaceAlpha', 0.6,'jitter','on','jitterAmount',0.15);
hold on;
scatter(5*ones(1,length(M2b_new.SP_CT)), M2b_new.SP_CT, dot_size, 'filled', 'MarkerFaceColor', col_blk, 'MarkerFaceAlpha', 0.6,'jitter','on','jitterAmount',0.15);
hold on;
scatter(7*ones(1,length(M3a_new.SP_CT)), M3a_new.SP_CT, dot_size, 'filled', 'MarkerFaceColor', col_blk, 'MarkerFaceAlpha', 0.6,'jitter','on','jitterAmount',0.15);
hold on;
scatter(8*ones(1,length(M3b_new.SP_CT)), M3b_new.SP_CT, dot_size, 'filled', 'MarkerFaceColor', col_blk, 'MarkerFaceAlpha', 0.6,'jitter','on','jitterAmount',0.15);
hold on;
scatter(10*ones(1,length(M5a_new.SP_CT)), M5a_new.SP_CT, dot_size, 'filled', 'MarkerFaceColor', col_blk, 'MarkerFaceAlpha', 0.6,'jitter','on','jitterAmount',0.15);
hold on;
scatter(11*ones(1,length(M5b_new.SP_CT)), M5b_new.SP_CT, dot_size, 'filled', 'MarkerFaceColor', col_blk, 'MarkerFaceAlpha', 0.6,'jitter','on','jitterAmount',0.15);
hold on;
plot([1, 4, 7, 10], meanM1235a, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 20);
hold on;
plot([2, 5, 8, 11], meanM1235b, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 20);
set(gca, 'Fontsize', 45);
box on;
xticks([1.5, 4.5, 7.5, 10.5]);
ylim([18, 45]);

set(findobj(gca,'type','line'),'linew',4)
set(gca,'linew',4)

xticks([1, 2, 4, 5, 7, 8, 10, 11]);
set(gca, 'XTickLabel', {'ABX  ', '  NoABX', 'ABX  ', '  NoABX', 'ABX  ', '  NoABX', 'ABX  ', '  NoABX'}, 'TickLabelInterpreter', 'latex', 'Fontsize', 38);

yticks([20, 25, 30, 35, 40, 45]);
set(gca, 'YTickLabel', {'20', ' 25', '30', ' 35', '40', '45'}, 'TickLabelInterpreter', 'latex', 'Fontsize', 45);

annotation('textbox', [0.14, 0.8, 0.1, 0.1], 'String', 'SP$\rightarrow$RSV', 'interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 45);

annotation('textbox', [0.33, 0.8, 0.1, 0.1], 'String', 'RSV \& SP', 'interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 45);

annotation('textbox', [0.53, 0.8, 0.1, 0.1], 'String', 'RSV$\rightarrow$SP', 'interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 45);

annotation('textbox', [0.73, 0.8, 0.1, 0.1], 'String', 'SP', 'interpreter', 'latex', 'EdgeColor', 'none', 'FontSize', 45);

xtickangle(30);

xline(3, 'LineWidth', 3);
xline(6, 'LineWidth', 3);
xline(9, 'LineWidth', 3);

yt = get(gca, 'YTick');
axis([xlim 0 ceil(max(yt) * 1.1)]);

xt = [10, 11];
yxis_reverse = 1;

if yxis_reverse == 0
    hold on;
    plot(xt([1 2]), [1 1] * max(yt) * 0.9, 'k-', 'LineWidth', 6);
    ylim([18, 45]);
end

dim = [.79 .59 .9 .3];
name = {'$ * $'};
annotation('textbox', dim, 'String', name, 'interpreter', 'latex', 'Fontsize', 55, 'Color', 'k', 'EdgeColor', 'none');

ylabel('SP Ct', 'Interpreter', 'latex', 'Fontsize', 60);

set(gca, 'Fontsize', 60);
box on;

if yxis_reverse == 1
    set(gca, 'YDir', 'reverse');
    ylim([15, 45]);
    hold on;
    plot(xt([1 2]), [1 1] * min(yt) * 0.96, 'k-', 'LineWidth', 6);
    hold on;
    yticks([20, 30, 40]);
    set(gca, 'YTickLabel', {'20', '30', '40'});
end


% saveas(gcf,'ABx_infants_v4.png')
