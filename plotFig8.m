% This script plots the first three colums in  fig 8

% clear all
clc
folderName = '/Users/meow/Documents/MATLAB/PH.D in UCI/multiple centroids_exp data & analysis/Data Analysis_mcmc/';
filename = [folderName,'subjData_mcmc_DeletingOutliers.mat'];
load(filename);
%%
numGroups = [3,3,4,4,4,4,5,6,6,7,8];
numItemsPerGroups = [10,10,8,6,6,6,6,5,5,4,4];

rmsError_multi = nan(3,11,2);
rmsError_singleton = nan(3,11,2);

meanError_multi = nan(3,11,2);
meanError_singleton = nan(3,11,2);

TW_multi = nan(3,11,2);
TW_singleton = nan(3,11,2);

filterRMS = nan(3,11,2);
filterMeanError = nan(3,11,2);


filterRMS_s = nan(3,11,2);
filterMeanError_s = nan(3,11,2);


cenrms = nan(3,11,2);
ratio = 0.05;

rmsMotorError = [14.8 14.6 13.4]';

for iExp = 1:11

   
    numColors = numGroups(iExp);
    N = numItemsPerGroups(iExp);
    totalItems = N*numColors;

    for iSub = 1:3
        subjectNum = iSub;
        
        if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_precue) 
            
            rmsError_multi(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).precue.errorRMS(:,1));
            rmsError_singleton(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).precue.errorRMS(:,2));

            meanError_multi(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).precue.meanerror(:,1));
            meanError_singleton(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).precue.meanerror(:,2));

            cenrms(iSub,iExp,1) = subjData_mcmc(iSub).experiment(iExp).precue.crms;
            TW_multi(iSub,iExp,1) = mean(diag(subjData_mcmc(iSub).experiment(iExp).precue.w));
            TW_singleton(iSub,iExp,1) = mean(diag(subjData_mcmc(iSub).experiment(iExp).precue.ws));
              
            filterRMS(iSub,iExp,1) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).precue.filterRMSError(:,1).^2));  
            filterMeanError(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).precue.filterMeanError(:,1));
            
            filterRMS_s(iSub,iExp,1) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).precue.filterRMSError(:,2).^2));  
            filterMeanError_s(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).precue.filterMeanError(:,2));

        end


         if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_postcue) 
            
            rmsError_multi(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.errorRMS(:,1));
            rmsError_singleton(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.errorRMS(:,2));

            meanError_multi(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,1));
            meanError_singleton(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,2));

            cenrms(iSub,iExp,2) = subjData_mcmc(iSub).experiment(iExp).precue.crms;
            TW_multi(iSub,iExp,2) = mean(diag(subjData_mcmc(iSub).experiment(iExp).postcue.w));
            TW_singleton(iSub,iExp,2) = mean(diag(subjData_mcmc(iSub).experiment(iExp).postcue.ws));
              
            filterRMS(iSub,iExp,2) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).postcue.filterRMSError(:,1).^2));
            filterMeanError(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.filterMeanError(:,1));
            
            filterRMS_s(iSub,iExp,2) = sqrt(mean(subjData_mcmc(iSub).experiment(iExp).postcue.filterRMSError(:,2).^2));
            filterMeanError_s(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).postcue.filterMeanError(:,2));

            

        end
         
      

    end

end

%%


expIdx = [2,3,7:10];

% post-cued multi filters
post_cued_f = TW_multi(:,expIdx,2);
newPost_cued_filter = [];
newPost_cued_filter = post_cued_f(:,1:3);
newPost_cued_filter(:,4) = mean(post_cued_f(:,4:5),2);
newPost_cued_filter(:,5) = post_cued_f(:,6);

% post-cued multi mean error magnitude
post_cued = meanError_multi(:,expIdx,2);
newPost_cued_error = [];
newPost_cued_error = post_cued(:,1:3);
newPost_cued_error(:,4) = mean(post_cued(:,4:5),2);
newPost_cued_error(:,5) = post_cued(:,6);

% post-cued singleton filters
post_cued_s  = meanError_singleton(:,expIdx,2);
newPost_cued_s_error = post_cued_s(:,1:3);
newPost_cued_s_error(:,4) = mean(post_cued_s(:,4:5),2);
newPost_cued_s_error(:,5) = post_cued_s(:,6);

% post-cued singleton mean error magnitude
post_cued_s_f = TW_singleton(:,expIdx,2);
newPost_cued_s_filter = post_cued_s_f(:,1:3);
newPost_cued_s_filter(:,4) = mean(post_cued_s_f(:,4:5),2);
newPost_cued_s_filter(:,5) = post_cued_s_f(:,6);



filter_meanError = [];
filter_meanError(:,1:3,:) = filterMeanError(:,expIdx(1:3),:);
filter_meanError(:,4,:) = mean(filterMeanError(:,expIdx(4:5),:),2);
filter_meanError(:,5,:) = filterMeanError(:,expIdx(end),:);


filter_meanError_s = [];
filter_meanError_s(:,1:3,:) = filterMeanError_s(:,expIdx(1:3),:);
filter_meanError_s(:,4,:) = mean(filterMeanError_s(:,expIdx(4:5),:),2);
filter_meanError_s(:,5,:) = filterMeanError_s(:,expIdx(end),:);



% pre-cued multi filters
Pre_cued_f = TW_multi(:,expIdx,1);
newPre_cued_filter = [];
newPre_cued_filter = Pre_cued_f(:,1:3);
newPre_cued_filter(:,4) = mean(Pre_cued_f(:,4:5),2);
newPre_cued_filter(:,5) = Pre_cued_f(:,6);

% pre-cued multi mean error magnitude
Pre_cued = meanError_multi(:,expIdx,1);
newPre_cued_error = [];
newPre_cued_error = Pre_cued(:,1:3);
newPre_cued_error(:,4) = mean(Pre_cued(:,4:5),2);
newPre_cued_error(:,5) = Pre_cued(:,6);

% pre-cued singleton filters
Pre_cued_s  = meanError_singleton(:,expIdx,1);
newPre_cued_s_error = [];
newPre_cued_s_error = Pre_cued_s(:,1:3);
newPre_cued_s_error(:,4) = mean(Pre_cued_s(:,4:5),2);
newPre_cued_s_error(:,5) = Pre_cued_s(:,6);

% pre-cued singleton mean error magnitude
Pre_cued_s_f = TW_singleton(:,expIdx,1);
newPre_cued_s_filter = [];
newPre_cued_s_filter = Pre_cued_s_f(:,1:3);
newPre_cued_s_filter(:,4) = mean(Pre_cued_s_f(:,4:5),2);
newPre_cued_s_filter(:,5) = Pre_cued_s_f(:,6);


excludeFlg  = 1; % subj 1: 3-7;, subj2&3: 3:5

if excludeFlg

   newPost_cued_filter(2:3,4:5) = nan;
   newPost_cued_s_filter(2:3,4:5) = nan;

   newPost_cued_error(2:3,4:5) = nan;
   newPost_cued_s_error(2:3,4:5) = nan;


   filter_meanError(2:3,4:5,2) = nan;
   filter_meanError_s(2:3,4:5,2) = nan;
end


%% fig 7 filters
color = [0.65 0.35 0.71;0.26 0.46 0.01];
markersize = 7;

figure(3)
subplot(2,5,1)

plot([1:5],mean(newPre_cued_filter,1),'o-','color',color(1,:),'LineWidth',1.5,'markersize',markersize);hold on

plot([1:5],mean(newPre_cued_s_filter,1),'o-','color',color(2,:),'LineWidth',1.5,'markersize',markersize);hold on

plot([1:5],ones(1,5),'-.k','LineWidth',1)

ylim([0.4,1.1])
xlim([0.5 5.5])


figure(3)
subplot(2,5,6)

plot([1:5],nanmean(newPost_cued_filter,1),'o-','color',color(1,:),'LineWidth',1.5,'markersize',markersize);hold on

plot([1:5],nanmean(newPost_cued_s_filter,1),'o-','color',color(2,:),'LineWidth',1.5,'markersize',markersize);hold on

plot([1:5],ones(1,5),'-.k','LineWidth',1)

ylim([0.4,1.1])
xlim([0.5 5.5])


%% figure 8 rms 

expIdx = [2,3,7,8,9,10];

multi_rms = [];
multi_rms(:,1:3,:) = rmsError_multi(:,expIdx(1:3),:);
multi_rms(:,4,:) = mean(rmsError_multi(:,expIdx(4:5),:),2);
multi_rms(:,5,:) = rmsError_multi(:,expIdx(end),:);

filter_rms = [];
filter_rms(:,1:3,:) = filterRMS(:,expIdx(1:3),:);
filter_rms(:,4,:) = sqrt(mean(filterRMS(:,expIdx(4:5),:).^2,2));
filter_rms(:,5,:) = filterRMS(:,expIdx(end),:);


filter_rms_s = [];
filter_rms_s(:,1:3,:) = filterRMS_s(:,expIdx(1:3),:);
filter_rms_s(:,4,:) = sqrt(mean(filterRMS_s(:,expIdx(4:5),:).^2,2));
filter_rms_s(:,5,:) = filterRMS_s(:,expIdx(end),:);


multi_meanE = [];
multi_meanE(:,1:3,:) = meanError_multi(:,expIdx(1:3),:);
multi_meanE(:,4,:) = mean(meanError_multi(:,expIdx(4:5),:),2);
multi_meanE(:,5,:) = meanError_multi(:,expIdx(end),:);


singleton_rms = [];
singleton_rms(:,1:3,:) = rmsError_singleton(:,expIdx(1:3),:);
singleton_rms(:,4,:) = mean(rmsError_singleton(:,expIdx(4:5),:),2);
singleton_rms(:,5,:) = rmsError_singleton(:,expIdx(end),:);

singleton_meanE = [];
singleton_meanE(:,1:3,:) = meanError_singleton(:,expIdx(1:3),:);
singleton_meanE(:,4,:) = mean(meanError_singleton(:,expIdx(4:5),:),2);
singleton_meanE(:,5,:) = meanError_singleton(:,expIdx(end),:);
 

cen_rms = [];

cen_rms(:,1:3,:) = cenrms(:,expIdx(1:3),:);
cen_rms(:,4,:) = mean(cenrms(:,expIdx(4:5),:),2);
cen_rms(:,5,:) = cenrms(:,expIdx(end),:);
 


t = length(expIdx)-1;

if excludeFlg

   multi_rms(2:3,4:5,2) = nan;
   singleton_rms(2:3,4:5,2) = nan;
   filter_rms(2:3,4:5,2) = nan;
   filter_rms_s(2:3,4:5,2) = nan;
end



figure(3)

subplot(2,5,2)

plot([1:5],nanmean(multi_rms(:,:,1).^2),'o-','color',color(1,:),'LineWidth',1.5,'markersize',markersize);hold on
plot([1:5],nanmean(singleton_rms(:,:,1).^2),'o-','color',color(2,:),'LineWidth',1.5,'markersize',markersize);hold on
plot([1:5],nanmean(cen_rms(:,:,1).^2),'k-.','LineWidth',1,'markersize',markersize); hold on

plot([1:5],nanmean(multi_rms(:,:,1).^2)-mean(rmsMotorError.^2),'o--','color',color(1,:),'LineWidth',1.5,'markersize',markersize);hold on
plot([1:5],nanmean(singleton_rms(:,:,1).^2)-mean(rmsMotorError.^2),'o--','color',color(2,:),'LineWidth',1.5,'markersize',markersize);hold on
plot([1:5],nanmean(cen_rms(:,:,2).^2),'k-.','LineWidth',1,'markersize',markersize); hold on




ylim([0, 31000])
xlim([0.5 5.5])
set(gca,'fontsize',14,'Ydir','reverse')


% ylim([0.4,1.1])
% xlim([0.5 5.5])
% 

figure(3)
subplot(2,5,7)

plot([1:5],nanmean(multi_rms(:,:,2).^2),'o-','color',color(1,:),'LineWidth',1.5,'markersize',markersize);hold on
plot([1:5],nanmean(singleton_rms(:,:,2).^2),'o-','color',color(2,:),'LineWidth',1.5,'markersize',markersize);hold on
plot([1:5],nanmean(cen_rms(:,:,2).^2),'k-.','LineWidth',1,'markersize',markersize); hold on

plot([1:5],nanmean(multi_rms(:,:,2).^2)-mean(rmsMotorError.^2),'o--','color',color(1,:),'LineWidth',1.5,'markersize',markersize);hold on
plot([1:5],nanmean(singleton_rms(:,:,2).^2)-mean(rmsMotorError.^2),'o--','color',color(2,:),'LineWidth',1.5,'markersize',markersize);hold on


ylim([0, 31000])
xlim([0.5 5.5])
set(gca,'fontsize',14,'Ydir','reverse')

%% fig 8 the third column

figure(3)
markersize = 10;


z = [nanmean(multi_rms(:,:,1).^2)'-nanmean(singleton_rms(:,:,1).^2)'-nanmean(filter_rms(:,:,1).^2)',nanmean(filter_rms(:,:,1).^2)',(nanmean(singleton_rms(:,:,1).^2)-mean(rmsMotorError.^2))',...
    ones(5,1)*mean(rmsMotorError.^2)];

z2 = [nanmean(multi_rms(:,:,2).^2)'-nanmean(singleton_rms(:,:,2).^2)'-nanmean(filter_rms(:,:,2).^2)',nanmean(filter_rms(:,:,2).^2)',(nanmean(singleton_rms(:,:,2).^2)-mean(rmsMotorError.^2))',...
    ones(5,1)*mean(rmsMotorError.^2)];

percentFlg = 1;
if percentFlg

    targetVariance = nanmean(cen_rms(:,:,2).^2);

    z_percent = (z./targetVariance');
    z2_percent = (z2./targetVariance');


     subplot(2,5,3)
    h = bar(1:5,z_percent,'stacked','BarWidth',0.6);
    h(1).FaceColor = [0.93 0.69 0.13];
    h(2).FaceColor = [0.85 0.33 0.1];
    h(3).FaceColor = [0.00,0.45,0.74];
    h(4).FaceColor = [0.65,0.65,0.65];


    h(1).EdgeColor = [0.93 0.69 0.13];
    h(2).EdgeColor = [0.85 0.33 0.1];
    h(3).EdgeColor = [0.00,0.45,0.74];
    h(4).EdgeColor = [0.65,0.65,0.65];

    xticks([1:5])
    xticklabels(num2str([3:7]'))
    legend('Channel+Centroid loss.','Filter loss','Working memory','Motor error')
    ylim([0 0.35])
    set(gca,'Ydir','reverse')
    set(gca,'FontSize',14)


    subplot(2,5,8)
    h = bar(1:5,z2_percent,'stacked','BarWidth',0.6);
    h(1).FaceColor = [0.93 0.69 0.13];
    h(2).FaceColor = [0.85 0.33 0.1];
    h(3).FaceColor = [0.00,0.45,0.74];
    h(4).FaceColor = [0.65,0.65,0.65];


    h(1).EdgeColor = [0.93 0.69 0.13];
    h(2).EdgeColor = [0.85 0.33 0.1];
    h(3).EdgeColor = [0.00,0.45,0.74];
    h(4).EdgeColor = [0.65,0.65,0.65];

    % xlim([0.5,5+0.5])
    xticks([1:5])
    xticklabels(num2str([3:7]'))
    set(gca,'Ydir','reverse')
    ylim([0 0.35])
    set(gca,'FontSize',14)



else

    subplot(2,5,3)
    h = bar(1:5,z,'stacked','BarWidth',0.6);
    h(1).FaceColor = [0.93 0.69 0.13];
    h(2).FaceColor = [0.85 0.33 0.1];
    h(3).FaceColor = [0.00,0.45,0.74];
    h(4).FaceColor = [0.65,0.65,0.65];

    h(1).EdgeColor = [0.93 0.69 0.13];
    h(2).EdgeColor = [0.85 0.33 0.1];
    h(3).EdgeColor = [0.00,0.45,0.74];
    h(4).EdgeColor = [0.65,0.65,0.65];
    xticks([1:5])
    xticklabels(num2str([3:7]'))
    legend('Perceptual loss','filter loss','Working memory','motor error')
    ylim([0, 18600])
    set(gca,'Ydir','reverse')
    set(gca,'FontSize',14)


    subplot(2,5,8)
    h = bar(1:5,z2,'stacked','BarWidth',0.6);
    h(1).FaceColor = [0.93 0.69 0.13];
    h(2).FaceColor = [0.85 0.33 0.1];
    h(3).FaceColor = [0.00,0.45,0.74];
    h(4).FaceColor = [0.65,0.65,0.65];

    h(1).EdgeColor = [0.93 0.69 0.13];
    h(2).EdgeColor = [0.85 0.33 0.1];
    h(3).EdgeColor = [0.00,0.45,0.74];
    h(4).EdgeColor = [0.65,0.65,0.65];

    % xlim([0.5,5+0.5])
    xticks([1:5])
    xticklabels(num2str([3:7]'))
    set(gca,'Ydir','reverse')
    ylim([0, 18600])
    set(gca,'FontSize',14)

end

plotItemAnalysis  % plot the last two colunms



