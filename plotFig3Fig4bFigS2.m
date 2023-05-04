% This script plots Fig 3, Fig 4b and Fig S2

clear
clc
close all

% load the data

filename = 'subjData_mcmc_DeletingOutliers.mat';
load(filename);
%%
numSubj = 3; 
numExp = 11;
numGroups = [3,3,4,4,4,4,5,6,6,7,8];

meanError =  nan(numSubj,numExp,4); % the third dimension: No distractors, precue, postcue, all-centroids 
targetFilter = nan(numSubj,numExp,4);
RMSError = nan(numSubj,numExp,4);


meanError_s =  nan(numSubj,numExp,4);
targetFilter_s = nan(numSubj,numExp,4);
RMSError_s = nan(numSubj,numExp,4);

Subj_efficiency = nan(numSubj,numExp,4);
Perfect_efficiency = nan(numSubj,numExp,4);


% regular trials
predictedMPlus1_TW = nan(numSubj,numExp,2);
predictedMPlus1_MeanError = nan(numSubj,numExp,2);

predictedMPlus2_TW = nan(numSubj,numExp);
predictedMPlus2_MeanError = nan(numSubj,numExp);

% singleton trials
predictedMPlus1_TW_S = nan(numSubj,numExp,2);
predictedMPlus1_MeanError_S = nan(numSubj,numExp,2);

predictedMPlus2_TW_S = nan(numSubj,numExp);
predictedMPlus2_MeanError_S = nan(numSubj,numExp);


ErrorSTD = nan(numSubj,numExp,4);
ErrorSTD_s = nan(numSubj,numExp,4); % singleton trials

SEforMotorError = nan(numSubj,numExp);
meanErrorforMotorError = nan(numSubj,numExp);

filterRMS = nan(numSubj,numExp,4);
filterRMS_s = nan(numSubj,numExp,4);
order = [1:3,7:11,4:6];
%%
for iExp = 1:11
    
    numColors = numGroups(iExp);
    
    for iSub = 1:3
        
        if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_Onetarget) % zero distactors
              
              % multi-items
              Subj_efficiency(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.efficiency)*numColors;
              Perfect_efficiency(iSub,iExp,1) = Subj_efficiency(iSub,iExp,1);
              targetFilter(iSub,iExp,1) =  mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.w(1,:));

              meanError(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.meanerror(:,1));
              RMSError(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.errorRMS(:,1));  
              ErrorSTD(iSub,iExp,1) = sqrt(sum(subjData_mcmc(iSub).experiment(iExp).Onetarget.StandardError(:,1).^2))/numColors;

              filterRMS(iSub,iExp,1) = subjData_mcmc(iSub).experiment(iExp).allCen.filterRMSErrorr ; 


              % singleton
              targetFilter_s(iSub,iExp,1) =  mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.ws(1,:));
              meanError_s(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.meanerror(:,2));
              RMSError_s(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.errorRMS(:,2));  
              ErrorSTD_s(iSub,iExp,1) = sqrt(sum(subjData_mcmc(iSub).experiment(iExp).Onetarget.StandardError(:,2).^2))/numColors;
              
        end
        
         if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_precue)  % precue


              Subj_efficiency(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).precue.subjEff.eff,'omitnan');
              Perfect_efficiency(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).precue.subjEff.effPerfect,'omitnan');
              targetFilter(iSub,iExp,2) = mean(diag(subjData_mcmc(iSub).experiment(iExp).precue.w));
              
              meanError(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).precue.meanerror(:,1));
              RMSError(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).precue.errorRMS(:,1));  
              ErrorSTD(iSub,iExp,2) = sqrt(sum(subjData_mcmc(iSub).experiment(iExp).precue.StandardError(:,1).^2))/numColors;
              
              filterRMS(iSub,iExp,2) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).precue.filterRMSError(:,1).^2));  


              targetFilter_s(iSub,iExp,2) = mean(diag(subjData_mcmc(iSub).experiment(iExp).precue.ws));              
              meanError_s(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).precue.meanerror(:,2));
              RMSError_s(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).precue.errorRMS(:,2));  
              ErrorSTD_s(iSub,iExp,2) = sqrt(sum(subjData_mcmc(iSub).experiment(iExp).precue.StandardError(:,2).^2))/numColors;
              
              filterRMS_s(iSub,iExp,2) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).precue.filterRMSError(:,2).^2));  

         end
         

        if isfield(subjData_mcmc(iSub).experiment(iExp),'allCen') && ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_Onetarget)

            Subj_efficiency(iSub,iExp,4) = subjData_mcmc(iSub).experiment(iExp).allCen.subjEff.eff;

            Perfect_efficiency(iSub,iExp,4) = subjData_mcmc(iSub).experiment(iExp).allCen.perfectEff.eff;

            meanError(iSub,iExp,4) = subjData_mcmc(iSub).experiment(iExp).allCen.meanerror(1);
            RMSError(iSub,iExp,4) = subjData_mcmc(iSub).experiment(iExp).allCen.errorRMS(1);
          
            targetFilter(iSub,iExp,4) = 1-var(subjData_mcmc(iSub).experiment(iExp).allCen.w);  % this is 1-var(filters)
            ErrorSTD(iSub,iExp,4) = subjData_mcmc(iSub).experiment(iExp).allCen.StandardError(1);
            
            filterRMS(iSub,iExp,4) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.filterRMSErrorr(:,1).^2));  
            filterRMS_s(iSub,iExp,4) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).Onetarget.filterRMSErrorr(:,2).^2));  

        end

         if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_postcue)

             Subj_efficiency(iSub,iExp,3) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.subjEff.eff,'omitnan');
             Perfect_efficiency(iSub,iExp,3) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.subjEff.effPerfect,'omitnan');

             targetFilter(iSub,iExp,3) = mean(diag(subjData_mcmc(iSub).experiment(iExp).postcue.w));
              meanError(iSub,iExp,3) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,1));
              RMSError(iSub,iExp,3) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.errorRMS(:,1));  
              ErrorSTD(iSub,iExp,3) = sqrt(sum(subjData_mcmc(iSub).experiment(iExp).postcue.StandardError(:,1).^2))/numColors;
              
              filterRMS(iSub,iExp,3) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).postcue.filterRMSError(:,1).^2));
                  
              targetFilter_s(iSub,iExp,3) = mean(diag(subjData_mcmc(iSub).experiment(iExp).postcue.ws));
              meanError_s(iSub,iExp,3) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,2));
              RMSError_s(iSub,iExp,3) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.errorRMS(:,2));  
              ErrorSTD_s(iSub,iExp,3) = sqrt(sum(subjData_mcmc(iSub).experiment(iExp).postcue.StandardError(:,2).^2))/numColors;
        
              filterRMS_s(iSub,iExp,3) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).postcue.filterRMSError(:,2).^2));
             


              % regular trials
              predictedMPlus1_TW(iSub,iExp,1) = subjData_mcmc(iSub).experiment(iExp).SummaryStats_postcue(3);
              predictedMPlus1_MeanError(iSub,iExp,1) =  subjData_mcmc(iSub).experiment(iExp).SummaryStats_postcue(1);
              
              
              predictedMPlus1_TW(iSub,iExp,2) = subjData_mcmc(iSub).experiment(iExp).PredictedMPlus1(1);
              predictedMPlus1_MeanError(iSub,iExp,2) =  subjData_mcmc(iSub).experiment(iExp).PredictedMPlus1(3);
              
              predictedMPlus2_TW(iSub,iExp) = subjData_mcmc(iSub).experiment(iExp).PredictedMPlus2(1);
              predictedMPlus2_MeanError(iSub,iExp) =  subjData_mcmc(iSub).experiment(iExp).PredictedMPlus2(3);
              
              % singleton trials
              
              predictedMPlus1_TW_S(iSub,iExp,1) = mean(diag(subjData_mcmc(iSub).experiment(iExp).postcue.ws));
              predictedMPlus1_MeanError_S(iSub,iExp,1) =  mean(subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,2));
              
              
              predictedMPlus1_TW_S(iSub,iExp,2) = subjData_mcmc(iSub).experiment(iExp).PredictedMPlus1Singleton(1);
              predictedMPlus1_MeanError_S(iSub,iExp,2) =  subjData_mcmc(iSub).experiment(iExp).PredictedMPlus1Singleton(3);
              
              predictedMPlus2_TW_S(iSub,iExp) = subjData_mcmc(iSub).experiment(iExp).PredictedMPlus2Singleton(1);
              predictedMPlus2_MeanError_S(iSub,iExp) =  subjData_mcmc(iSub).experiment(iExp).PredictedMPlus2Singleton(3);
              
        end
         
                
    end
    
end


%% mean error magnitude

figure(1)
logflg = 1;

subplot(3,2,1)
interval = [-0.1,0,0.1];
marker_avg = {'-','-','-'};
marker_subj = {'d','o','<'};
color = [1,0.3250,0.0980;0.9290,0.6940,0.1250;0.3010,0.7450,0.9330;0.45,0.99,0.07];
% color = [1,0.3250,0.0980;0.9290,0.6940,0.1250;0.45,0.99,0.07];


if logflg
    meanError_newOrder = log10(meanError(:,order,:));
else
    meanError_newOrder = meanError(:,order,:);
end

failedIndex{1} = 8;
failedIndex{2} = [5,6,7];
failedIndex{3} = [5,6,7];


markerSize1 = 8;
markerSize2 = 5;

xindex = [1:11];
linewidth = 1.5;

for isub = 1:3
    
    xindexSc = xindex;
    xindexSc(failedIndex{isub})=[];
     
    plot(xindex+interval(isub),meanError_newOrder(isub,:,1),marker_subj{isub},'color',color(1,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    plot(xindex+interval(isub),meanError_newOrder(isub,:,2),marker_subj{isub},'color',color(2,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    
    plot(xindexSc+interval(isub),meanError_newOrder(isub,xindexSc,3),marker_subj{isub},'color',color(3,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    plot(failedIndex{isub}+interval(isub),meanError_newOrder(isub,failedIndex{isub},3),marker_subj{isub},'color',color(3,:),'markersize',markerSize2 ,'LineWidth',linewidth);hold on

    plot(xindex+interval(isub),meanError_newOrder(isub,:,4),marker_subj{isub},'color',color(4,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
     
       
end

g1 = 1:8;
g2 = 9:11;


plot(g1,mean(meanError_newOrder(:,g1,1),'omitnan'),marker_avg{1},'color',color(1,:),'LineWidth',linewidth);hold on
plot(g1,mean(meanError_newOrder(:,g1,2),'omitnan'),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g1,mean(meanError_newOrder(:,g1,3),'omitnan'),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on
plot(g1,mean(meanError_newOrder(:,g1,4),'omitnan'),marker_avg{3},'color',color(4,:),'LineWidth',linewidth);hold on


plot(g2,mean(meanError_newOrder(:,g2,1)),marker_avg{1},'color',color(1,:),'LineWidth',linewidth);hold on
plot(g2,mean(meanError_newOrder(:,g2,2)),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g2,mean(meanError_newOrder(:,g2,3)),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on
plot(g2,mean(meanError_newOrder(:,g2,4),'omitnan'),marker_avg{3},'color',color(4,:),'LineWidth',linewidth);hold on

xticks(1:11)
xlim([0,11+1])
% t = {'C3S0D10 NR \newline    3';'C3S0D10 \newline    3';'C4S0D8 \newline    4';'C4S0D6 \newline    4';'C2S2D6 \newline    4';'C0S4D6 \newline    4';...
%     'C5S0D6 \newline    5' ;'C6S0D5 \newline    6 ';'C3S3D5 \newline    6';'C4S3D4 \newline    7';'C4S4D4 \newline    8'};
% t = t(order);
xticklabels([])
ylabel('Mean error magnitude')
set(gca,'fontsize',16,'FontName','Arial')
set(gca,'LineWidth',linewidth)

% legend('precue','postcue','precue filters after decimation','postcue filter error after decimation')

% legend('Subj1 No-distractors','Subj1 precue','Subj1 postcue','Subj1 ALL','Subj2 No distractors','Subj2 precue','Subj2 postcue','Subj2 ALL',...
%     'Subj3 No distractors','Subj3 precue','Subj3 postcue','Subj3 ALL','mean zero-distractor','mean precue','mean postcue','mean ALL')

stdError = 2*squeeze(sqrt(nansum(ErrorSTD(:,order,:).^2))./sum(~isnan(ErrorSTD(:,order,:))));

std_g1 = mean(stdError(g1,1));
if logflg ==1
    tt_g1  = mean(10.^meanError_newOrder(:,g1,1),'omitnan');
    std_g1 = log10(mean(tt_g1)+std_g1)-log10(mean(tt_g1));
     errorbar(g1(end)+0.3,log10(mean(tt_g1)),std_g1,marker_avg{1},'color',color(1,:),'LineWidth',linewidth,'marker','.','markersize',10);hold on

else
    tt_g1  = mean(meanError_newOrder(:,g1,1),'omitnan');
    errorbar(g1(end)+0,3,mean(tt_g1),std_g1,marker_avg{1},'color',color(1,:),'LineWidth',linewidth,'marker','.','markersize',10);hold on

end
 
std_g1 = mean(stdError(g1,4));
if logflg ==1
    tt_g1  = mean(10.^meanError_newOrder(:,g1,4),'omitnan');
    std_g1 = log10(mean(tt_g1)+std_g1)-log10(mean(tt_g1));
    errorbar(g1(end)+0.5,log10(mean(tt_g1)),std_g1,marker_avg{1},'color',color(4,:),'LineWidth',linewidth,'marker','.','markersize',10);hold on

else
    tt_g1  = mean(meanError_newOrder(:,g1,4),'omitnan');
  errorbar(g1(end)+0.5,mean(tt_g1),std_g1,marker_avg{1},'color',color(4,:),'LineWidth',linewidth,'marker','.','markersize',10);hold on

end


std_g2 = mean(stdError(g2,1));
if logflg ==1
    tt_g2  = mean(10.^meanError_newOrder(:,g2,1),'omitnan');
    std_g2 = log10(mean(tt_g2)+std_g2)-log10(mean(tt_g2));
     errorbar(g2(end)+0.3,log10(mean(tt_g2)),std_g2,marker_avg{1},'color',color(1,:),'LineWidth',linewidth,'marker','.','markersize',10);hold on

else
    tt_g2  = mean(meanError_newOrder(:,g2,1),'omitnan');
    errorbar(g2(end)+0.3,mean(tt_g2),std_g2,marker_avg{1},'color',color(1,:),'LineWidth',linewidth,'marker','.','markersize',10);hold on

end
 
std_g2 = mean(stdError(g2,4));
if logflg ==1
    tt_g2  = mean(10.^meanError_newOrder(:,g2,4),'omitnan');
    std_g2 = log10(mean(tt_g2)+std_g2)-log10(mean(tt_g2));
    errorbar(g2(end)+0.5,log10(mean(tt_g2)),std_g2,marker_avg{1},'color',color(4,:),'LineWidth',linewidth,'marker','.','markersize',10);hold on

else
    tt_g2  = mean(meanError_newOrder(:,g2,4),'omitnan');
  errorbar(g2(end)+0.5,mean(tt_g2),std_g2,marker_avg{1},'color',color(4,:),'LineWidth',linewidth,'marker','.','markersize',10);hold on

end
% 
edgealpha = 0.5;
% 
if logflg

    y = log10([mean(10.^meanError_newOrder(:,g1,2),'omitnan')-stdError(g1,2)', fliplr(mean(10.^meanError_newOrder(:,g1,2),'omitnan')+stdError(g1,2)')]);
else
    y = [mean(meanError_newOrder(:,g1,2),'omitnan')-stdError(g1,2)', fliplr(mean(meanError_newOrder(:,g1,2),'omitnan')+stdError(g1,2)')];
end
h = fill([g1,fliplr(g1)],y,color(2,:));hold on
h.EdgeAlpha = edgealpha ;
h.FaceAlpha = 0.2;
h.EdgeColor = color(2,:);   

if logflg

    y = log10([mean(10.^meanError_newOrder(:,g1,3),'omitnan')-stdError(g1,3)', fliplr(mean(10.^meanError_newOrder(:,g1,3),'omitnan')+stdError(g1,3)')]);
else
    y = [mean(meanError_newOrder(:,g1,3),'omitnan')-stdError(g1,3)', fliplr(mean(meanError_newOrder(:,g1,3),'omitnan')+stdError(g1,3)')];
end

h = fill([g1,fliplr(g1)],y,color(3,:));hold on
h.EdgeAlpha = edgealpha;
h.FaceAlpha = 0.2;
h.EdgeColor = color(3,:);

      
if logflg
    y = log10([mean(10.^meanError_newOrder(:,g2,2),'omitnan')-stdError(g2,2)', fliplr(mean(10.^meanError_newOrder(:,g2,2),'omitnan')+stdError(g2,2)')]);
else
    y = [mean(meanError_newOrder(:,g2,2),'omitnan')-stdError(g2,2)', fliplr(mean(meanError_newOrder(:,g2,2),'omitnan')+stdError(g2,2)')];
end
h = fill([g2,fliplr(g2)],y,color(2,:));hold on

h.EdgeAlpha = edgealpha ;
h.FaceAlpha = 0.2;
h.EdgeColor = color(2,:);


if logflg
    y = log10([mean(10.^meanError_newOrder(:,g2,3),'omitnan')-stdError(g2,3)', fliplr(mean(10.^meanError_newOrder(:,g2,3),'omitnan')+stdError(g2,3)')]);
else
    y = [mean(meanError_newOrder(:,g2,3),'omitnan')-stdError(g2,3)', fliplr(mean(meanError_newOrder(:,g2,3),'omitnan')+stdError(g2,3)')];
end
h = fill([g2,fliplr(g2)],y,color(3,:));hold on
h.EdgeAlpha = edgealpha;
h.FaceAlpha = 0.2;
h.EdgeColor = color(3,:);

legend off

if logflg
    t = [10,20,30,40,80,120];
    yticks(log10(t))
    ylim([0.75 2.2])
    yticklabels(t)
else
    ylim([0 120])
end

set(gca,'Ydir','reverse')

%% target filters

g1 = 1:8;
g2 = 9:11;


makersize1 = 7;
figure(1);
% h.Position = [12.00 291.00 1375.00 575.00];
subplot(3,2,3)
interval = [-0.1,0,0.1];
marker_avg = {'-','-','-','-'};
marker_subj = {'d','o','<'};
color = [1,0.3250,0.0980;0.9290,0.6940,0.1250;0.3010,0.7450,0.9330;0.45,0.99,0.07];
% color = [0.3010,0.7450,0.9330;1,0.3250,0.0980;0.9290,0.6940,0.1250;0.11,0.41,0.11];
targetFilter(1,:,1) = 1;
targetFilter(2:3,1:10,1) = 1;
targetFilter_newOrder = targetFilter(:,order,:);

for isub = 1:3

    xindexSc = xindex;
    xindexSc(failedIndex{isub})=[];
   

    plot([1:11]+interval(isub),targetFilter_newOrder(isub,:,1),marker_subj{isub},'color',color(1,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    plot([1:11]+interval(isub),targetFilter_newOrder(isub,:,2),marker_subj{isub},'color',color(2,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    
    plot(xindexSc+interval(isub),targetFilter_newOrder(isub,xindexSc,3),marker_subj{isub},'color',color(3,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    plot(failedIndex{isub}+interval(isub),targetFilter_newOrder(isub,failedIndex{isub},3),marker_subj{isub},'color',color(3,:),'markersize',markerSize2 ,'LineWidth',linewidth);hold on

      
    plot([1:11]+interval(isub),targetFilter_newOrder(isub,:,4),marker_subj{isub},'color',color(4,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on   
 
end
%

plot(g1,mean(targetFilter_newOrder(:,g1,1),'omitnan'),marker_avg{1},'color',color(1,:),'LineWidth',linewidth);hold on
plot(g1,mean(targetFilter_newOrder(:,g1,2),'omitnan'),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g1,mean(targetFilter_newOrder(:,g1,3),'omitnan'),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on
 plot(g1,mean(targetFilter_newOrder(:,g1,4),'omitnan'),marker_avg{4},'color',color(4,:),'LineWidth',linewidth);hold on

plot(g2,mean(targetFilter_newOrder(:,g2,1),'omitnan'),marker_avg{1},'color',color(1,:),'LineWidth',linewidth);hold on
plot(g2,mean(targetFilter_newOrder(:,g2,2),'omitnan'),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g2,mean(targetFilter_newOrder(:,g2,3),'omitnan'),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on
plot(g2,mean(targetFilter_newOrder(:,g2,4),'omitnan'),marker_avg{4},'color',color(4,:),'LineWidth',linewidth);hold on

ylim([0.2,1.1])
yticks([0.2:0.2:1.1]);
ylabel('TW')
set(gca, 'YColor', 'k');
set(gca,'fontsize',16,'FontName','Arial')


xlim([0,11+1])
xticks([1:11])
yyaxis right
ylim([0.2,1.1])
ytick = [0.2:0.2:1.1];
yticks(ytick)
t = ytick./(1-ytick);
yticklabels([num2str(round(t',2))])
ylabel('TW/(1-TW)')
set(gca, 'YColor', 'k');
set(gca,'fontsize',16,'FontName','Arial')
set(gca,'LineWidth',linewidth)
% title('target weight')
xticklabels([])

%% singleton mean error magnitude
% close all
figure(1) % mean error magnitude

subplot(3,2,2)
color = [1,0.3250,0.0980;0.9290,0.6940,0.1250;0.3010,0.7450,0.9330;0.45,0.99,0.07];
marker_subj = {'d','o','<'};
marker_avg = {'-','-','-'};
interval = [-0.1,0,0.1];
meanError_s_newOrder = log10(meanError_s(:,order,:));
for isub = 1:3

    xindexSc = xindex;
    xindexSc(failedIndex{isub})=[];
  
    
    plot([1:11]+interval(isub),meanError_s_newOrder(isub,:,1),marker_subj{isub},'color',color(1,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    plot([1:11]+interval(isub),meanError_s_newOrder(isub,:,2),marker_subj{isub},'color',color(2,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    
    plot(xindexSc+interval(isub),meanError_s_newOrder(isub,xindexSc,3),marker_subj{isub},'color',color(3,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    plot(failedIndex{isub}+interval(isub),meanError_s_newOrder(isub,failedIndex{isub},3),marker_subj{isub},'color',color(3,:),'markersize',markerSize2 ,'LineWidth',linewidth);hold on

%     plot([1:11]+interval(isub),meanError_s_newOrder(isub,:,3),marker_subj{isub},'color',color(3,:),'markersize',makersize1 ,'LineWidth',linewidth);hold on

end
 
plot(g1,nanmean(meanError_s_newOrder(:,g1,1)),marker_avg{1},'color',color(1,:),'LineWidth',linewidth);hold on
plot(g1,nanmean(meanError_s_newOrder(:,g1,2)),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g1,nanmean(meanError_s_newOrder(:,g1,3)),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on


plot(g2,nanmean(meanError_s_newOrder(:,g2,1)),marker_avg{1},'color',color(1,:),'LineWidth',linewidth);hold on
plot(g2,nanmean(meanError_s_newOrder(:,g2,2)),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g2,nanmean(meanError_s_newOrder(:,g2,3)),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on


stdError = 2*squeeze(sqrt(nansum(ErrorSTD_s(:,order,:).^2))./sum(~isnan(ErrorSTD_s(:,order,:))));

std_g1 = mean(stdError(g1,1));
if logflg ==1
    tt_g1  = mean(10.^meanError_s_newOrder(:,g1,1),'omitnan');
    std_g1 = log10(mean(tt_g1)+std_g1)-log10(mean(tt_g1));
    errorbar(g1(end)+0.3,log10(mean(tt_g1)),std_g1,marker_avg{1},'color',color(1,:),'LineWidth',linewidth,'marker','.','markersize',makersize1 );hold on

else
    tt_g1  = mean(meanError_s_newOrder(:,g1,1),'omitnan');
    errorbar(g1(end)+0.3,mean(tt_g1),std_g1,marker_avg{1},'color',color(1,:),'LineWidth',linewidth,'marker','.','markersize',makersize1 );hold on

end
 

std_g2 = mean(stdError(g2,1));
if logflg ==1
    tt_g2  = mean(10.^meanError_s_newOrder(:,g2,1),'omitnan');
    std_g2 = log10(mean(tt_g2)+std_g2)-log10(mean(tt_g2));
     errorbar(g2(end)+0.3,log10(mean(tt_g2)),std_g2,marker_avg{1},'color',color(1,:),'LineWidth',linewidth,'marker','.','markersize',makersize1 );hold on

else
    tt_g2  = mean(meanError_s_newOrder(:,g2,1),'omitnan');
    errorbar(g2(end)+0.3,mean(tt_g2),std_g2,marker_avg{1},'color',color(1,:),'LineWidth',linewidth,'marker','.','markersize',makersize1 );hold on

end

edgealpha = 0.5;

if logflg

    y = log10([mean(10.^meanError_s_newOrder(:,g1,2),'omitnan')-stdError(g1,2)', fliplr(mean(10.^meanError_s_newOrder(:,g1,2),'omitnan')+stdError(g1,2)')]);
else
    y = [mean(meanError_s_newOrder(:,g1,2),'omitnan')-stdError(g1,2)', fliplr(mean(meanError_s_newOrder(:,g1,2),'omitnan')+stdError(g1,2)')];
end
h = fill([g1,fliplr(g1)],y,color(2,:));hold on
h.EdgeAlpha = edgealpha ;
h.FaceAlpha = 0.2;
h.EdgeColor = color(2,:);   

if logflg

    y = log10([mean(10.^meanError_s_newOrder(:,g1,3),'omitnan')-stdError(g1,3)', fliplr(mean(10.^meanError_s_newOrder(:,g1,3),'omitnan')+stdError(g1,3)')]);
else
    y = [mean(meanError_s_newOrder(:,g1,3),'omitnan')-stdError(g1,3)', fliplr(mean(meanError_s_newOrder(:,g1,3),'omitnan')+stdError(g1,3)')];
end

h = fill([g1,fliplr(g1)],y,color(3,:));hold on
h.EdgeAlpha = edgealpha;
h.FaceAlpha = 0.2;
h.EdgeColor = color(3,:);

if logflg
    y = log10([mean(10.^meanError_s_newOrder(:,g2,2),'omitnan')-stdError(g2,2)', fliplr(mean(10.^meanError_s_newOrder(:,g2,2),'omitnan')+stdError(g2,2)')]);
else
    y = [mean(meanError_s_newOrder(:,g2,2),'omitnan')-stdError(g2,2)', fliplr(mean(meanError_s_newOrder(:,g2,2),'omitnan')+stdError(g2,2)')];
end
h = fill([g2,fliplr(g2)],y,color(2,:));hold on

h.EdgeAlpha = edgealpha ;
h.FaceAlpha = 0.2;
h.EdgeColor = color(2,:);


if logflg
    y = log10([mean(10.^meanError_s_newOrder(:,g2,3),'omitnan')-stdError(g2,3)', fliplr(mean(10.^meanError_s_newOrder(:,g2,3),'omitnan')+stdError(g2,3)')]);
else
    y = [mean(meanError_s_newOrder(:,g2,3),'omitnan')-stdError(g2,3)', fliplr(mean(meanError_s_newOrder(:,g2,3),'omitnan')+stdError(g2,3)')];
end
h = fill([g2,fliplr(g2)],y,color(3,:));hold on
h.EdgeAlpha = edgealpha;
h.FaceAlpha = 0.2;
h.EdgeColor = color(3,:);


legend off

if logflg
    t = [10,20,30,40,80,120];
    yticks(log10(t))
    ylim([0.75 2.2])
    yticklabels(t)
%     title('multi-item log10 scale')
else
    ylim([0 120])
%     title('multi-item linear scale')
end


xticks(1:11)
xlim([0,11+1])
% t = {'C3S0D10 NR \newline    3';'C3S0D10 \newline    3';'C4S0D8 \newline    4';'C4S0D6 \newline    4';'C2S2D6 \newline    4';'C0S4D6 \newline    4';...
%     'C5S0D6 \newline    5' ;'C6S0D5 \newline    6 ';'C3S3D5 \newline    6';'C4S3D4 \newline    7';'C4S4D4 \newline    8'};
% t = t(order);
xticklabels([])
% title ('Singleton trials: Mean Error Magnitude')
ylabel('mean error magnitude (pixels)')
set(gca,'fontsize',16,'FontName','Arial')
set(gca,'LineWidth',linewidth)

legend('Subj1 No-distractor','Subj1 pre-cued','Subj1 post-cued','Subj2 No-distractor','Subj2 pre-cued','Subj2 post-cued','Subj1 No-distractor','Subj3 pre-cued','Subj3 post-cued','No-distractor average','pre-cued average','post-cued average','Location','best')


legend off
t = [10,20,40,80,120];
yticks(log10(t))
ylim([0.75 2.2])

yticklabels(t)
set(gca,'Ydir','reverse')

%% singleton target filters
targetWeight_s_newOrder = targetFilter_s(:,order,:);
figure(1)
subplot(3,2,4)
color = [1,0.3250,0.0980;0.9290,0.6940,0.1250;0.3010,0.7450,0.9330;0.45,0.99,0.07];
for isub = 1:3
    
    xindexSc = xindex;
    xindexSc(failedIndex{isub})=[];
  

    plot([1:11]+interval(isub),targetWeight_s_newOrder(isub,:,2),marker_subj{isub},'color',color(2,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
%     plot([1:11]+interval(isub),targetWeight_s_newOrder(isub,:,3),marker_subj{isub},'color',color(3,:),'markersize',makersize1 ,'LineWidth',linewidth);hold on
    
    plot(xindexSc+interval(isub),targetWeight_s_newOrder(isub,xindexSc,3),marker_subj{isub},'color',color(3,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    plot(failedIndex{isub}+interval(isub),targetWeight_s_newOrder(isub,failedIndex{isub},3),marker_subj{isub},'color',color(3,:),'markersize',markerSize2 ,'LineWidth',linewidth);hold on

    
end

plot(g1,mean(targetWeight_s_newOrder(:,g1,2)),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g1,mean(targetWeight_s_newOrder(:,g1,3)),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on

plot(g2,mean(targetWeight_s_newOrder(:,g2,2)),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g2,mean(targetWeight_s_newOrder(:,g2,3)),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on

legend('Subj1 pre-cued','Subj1 post-cued','Subj2 pre-cued','Subj2 post-cued','Subj3 pre-cued','Subj3 post-cued','pre-cued average','post-cued average','Location','best')


xticks(1:11)
xlim([0,11+1])
xticklabels([])
% xticklabels({'C3S0D10 NR \newline    3';'C3S0D10 \newline    3';'C4S0D8 \newline    4';'C4S0D6 \newline    4';'C2S2D6 \newline    4';'C0S4D6 \newline    4';...
%     'C5S0D6 \newline    5' ;'C6S0D5 \newline    6 ';'C3S3D5 \newline    6';'C4S3D4 \newline    7';'C4S4D4 \newline    8'})
% 
% title ('singleton trials')
ylim([0.2,1.1])
ytick = 0.2:0.2:1.1;
yticks(ytick)
% title('Singleton trials: Target Weight (TW)')
ylabel('Target Weight (TW)')
set(gca,'fontsize',16,'FontName','Arial')
set(gca,'LineWidth',linewidth)

yyaxis right
ylim([0.2,1.1])
yticks(ytick)
t = ytick./(1-ytick);
yticklabels([num2str(round(t',2))])
ylabel('TW/(1-TW)')
set(gca, 'YColor', 'k');
set(gca,'fontsize',16,'FontName','Arial')
set(gca,'LineWidth',linewidth)
legend off

%%
figure(1);
% h.Position = [34.00 376.00 1393.00 489.00];
subplot(3,2,5)

color = [1,0.3250,0.0980;0.9290,0.6940,0.1250;0.3010,0.7450,0.9330;0.45,0.99,0.07];
marker_subj = {'d','o','<'};
marker_avg = {'-','-','-'};
interval = [-0.1,0,0.1];

Subj_efficiency_newOrder =  Subj_efficiency(:,order,:);

Perfect_efficiency_newOrder = Perfect_efficiency(:,order,:);

for isub = 1:3
    

    xindexSc = xindex;
    xindexSc(failedIndex{isub})=[];
  
    plot([1:11]+interval(isub),Subj_efficiency_newOrder(isub,:,1).*numGroups,marker_subj{isub},'color',color(1,:),'markersize',markerSize1,'LineWidth',linewidth);hold on
    plot([1:11]+interval(isub),Subj_efficiency_newOrder(isub,:,2),marker_subj{isub},'color',color(2,:),'markersize',markerSize1,'LineWidth',linewidth );hold on
    

    plot(xindexSc+interval(isub),Subj_efficiency_newOrder(isub,xindexSc,3),marker_subj{isub},'color',color(3,:),'markersize',markerSize1 ,'LineWidth',linewidth);hold on
    plot(failedIndex{isub}+interval(isub),Subj_efficiency_newOrder(isub,failedIndex{isub},3),marker_subj{isub},'color',color(3,:),'markersize',markerSize2 ,'LineWidth',linewidth);hold on
  
    plot([1:11]+interval(isub),Subj_efficiency_newOrder(isub,:,4),marker_subj{isub},'color',color(4,:),'markersize',markerSize1,'LineWidth',linewidth);hold on


end
 

xticks(1:11)
xlim([0,11+1])
% xticklabels({'C3S0D10 NR \newline    3';'C3S0D10 \newline    3';'C4S0D8 \newline    4';'C4S0D6 \newline    4';'C2S2D6 \newline    4';'C0S4D6 \newline    4';...
%     'C5S0D6 \newline    5' ;'C6S0D5 \newline    6 ';'C3S3D5 \newline    6';'C4S3D4 \newline    7';'C4S4D4 \newline    8'})
xticklabels([])
set(gca,'fontsize',16,'FontName','Arial')
set(gca,'LineWidth',linewidth)
ylim([0,32])
ylabel('Subject''s Efficiency ')
% ylim([0,32])
% ylabel('Ideal detector Efficiency ') '

% 
plot(g1,mean(Subj_efficiency_newOrder(:,g1,1),'omitnan'),marker_avg{1},'color',color(1,:),'LineWidth',linewidth);hold on
plot(g1,mean(Subj_efficiency_newOrder(:,g1,2),'omitnan'),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g1,mean(Subj_efficiency_newOrder(:,g1,3),'omitnan'),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on
plot(g1,mean(Subj_efficiency_newOrder(:,g1,4),'omitnan'),marker_avg{3},'color',color(4,:),'LineWidth',linewidth);hold on



plot(g2,mean(Subj_efficiency_newOrder(:,g2,1),'omitnan'),marker_avg{1},'color',color(1,:),'LineWidth',linewidth);hold on
plot(g2,mean(Subj_efficiency_newOrder(:,g2,2),'omitnan'),marker_avg{2},'color',color(2,:),'LineWidth',linewidth);hold on
plot(g2,mean(Subj_efficiency_newOrder(:,g2,3),'omitnan'),marker_avg{3},'color',color(3,:),'LineWidth',linewidth);hold on
plot(g2,mean(Subj_efficiency_newOrder(:,g2,4),'omitnan'),marker_avg{3},'color',color(4,:),'LineWidth',linewidth);hold on

% 
% 
effiency_averaged_postcue = nanmean(Perfect_efficiency_newOrder(:,:,3));
plot(g1,effiency_averaged_postcue(g1),'.-','color',[0,0,1],'markersize',makersize1 ,'LineWidth',linewidth); hold on
plot(g2,effiency_averaged_postcue(g2),'.-','color',[0,0,1],'markersize',makersize1 ,'LineWidth',linewidth); hold on

legend off

%% plot figure 4b / S2

fig4Flg = 1; % 1 plot figure 4b; else plot figure S2

Perfect_efficiency_successed = Perfect_efficiency;
Perfect_efficiency_successed(2:3,8:10,3) = nan;
Perfect_efficiency_successed(1,11,3) = nan;

Perfect_efficiency_failed = nan(size(Perfect_efficiency));
Perfect_efficiency_failed(2:3,8:10,3) = Perfect_efficiency(2:3,8:10,3);
Perfect_efficiency_failed(1,11,3) = Perfect_efficiency(1,11,3);



meanError_successed = meanError;
meanError_successed(2:3,8:10,3) = nan;
meanError_successed(1,11,3) = nan;

meanError_failed = nan(size(meanError));
meanError_failed(2:3,8:10,3) = meanError(2:3,8:10,3);
meanError_failed(1,11,3) = meanError(1,11,3);


meanError_subj = squeeze(nanmean(meanError_successed(:,:,1:3),1));
meanError_fail = squeeze(nanmean(meanError_failed(:,:,1:3),1));

NumOfTargetProcessedSuccessed = squeeze(nanmean(Perfect_efficiency_successed(:,:,1:3),1)./numGroups);

numItemsPerGroups = [10,10,8,6,6,6,6,5,5,4,4];
expfolder = '/Users/meow/Documents/MATLAB/PH.D in UCI/multiple centroids_exp data & analysis/exp*';
expfolders = dir(expfolder);
p = 0.01:0.01:0.99;
h = figure();
ExpIdx = [1,2,3,4,7,8, 10,11];
colors = [0.97,0.61,0.87;0.39,0.83,0.07;0.99,0.52,0.23;0.00,0.00,1.00;1.00,0.07,0.65;0.31,0.45,0.14;0.06,1.00,1.00; 1,0,0]';
z = [0.03 0.56 0.88;0.67 0.19 0.37;0.46, 0.98 0.16]';
colors2(:,1:4) = colors(:,1:4);
colors2(:,5:6) =  z(:,1:2);
colors2(:,7) = colors(:,5);
colors2(:,8) = z(:,3);
colors2(:,9:11) = colors(:,6:8);


for i = 1:8
    iExp = ExpIdx(i);

    selectedfolder = [expfolders(iExp).folder,'/',expfolders(iExp).name];
    numColors = numGroups(iExp);
    N = numItemsPerGroups(iExp);

    if fig4Flg % plot Fig 4b
        dots = N-N*p;
    else   % plot Fig S2

        items = N*numColors;
        dots = items-items*p;
    end



    efficiencyFile = ['EFF_fxdD100_RovStd80*.mat'];
    file = dir([selectedfolder,'/',efficiencyFile]);
    load([file.folder,'/',file.name]);
    
    if i  == 4 || i==8
     plot(data(:,1),dots,'color',colors(:,i),'linestyle','-','LineWidth',linewidth); hold on
    else
     plot(data(:,1),dots,'color',colors(:,i),'linestyle','-','LineWidth',linewidth); hold on
    end

end

eff_successed = nan(11,3);
eff_failed = nan(11,3);

for iExp = 1:11

    selectedfolder = [expfolders(iExp).folder,'/',expfolders(iExp).name];
    numColors = numGroups(iExp);
    N = numItemsPerGroups(iExp);
    if fig4Flg % plot Fig 4b
        dots = N-N*p;
    else   % plot Fig S2

        items = N*numColors;
        dots = items-items*p;
    end
    efficiencyFile = ['EFF_fxdD100_RovStd80*.mat'];
    file = dir([selectedfolder,'/',efficiencyFile]);
    load([file.folder,'/',file.name]);
     
    for i = 1:3
        if ~isnan(meanError_subj(iExp,i))

            [~,idx_precue] = min(abs(data(:,1)-meanError_subj(iExp,i)));
            eff_successed(iExp,i) = dots(idx_precue);
        else
            eff_successed(iExp,i) = nan;
        end
    end


    for i = 1:3
        if ~isnan(meanError_fail(iExp,i))

            [~,idx_precue] = min(abs(data(:,1)-meanError_fail(iExp,i)));
            eff_failed(iExp,i) = dots(idx_precue);
        else
            eff_failed(iExp,i) = nan;
        end
    end

                   
end

for i = 1:11
  plot(meanError_subj(i,1),eff_successed(i,1),"+",'color',colors2(:,i),'MarkerSize',10,'linewidth',2); hold on
  plot(meanError_subj(i,2),eff_successed(i,2),"x",'color',colors2(:,i),'MarkerSize',10,'linewidth',2); hold on

  if ~isnan(meanError_subj(i,3))
      plot(meanError_subj(i,3),eff_successed(i,3),"*",'color',colors2(:,i),'MarkerSize',10,'linewidth',2); hold on
  end
  
  if ~isnan(meanError_fail(i,3))
      plot(meanError_fail(i,3),eff_failed(i,3),'square','color',colors2(:,i),'MarkerSize',10,'linewidth',2); hold on
  end
end

xlabel('mean error magnitude')

if fig4Flg % plot Fig 4b
    ylabel('the number of surviving target items')
else   % plot Fig S2

    ylabel('the number of stimulus items processed')
end



