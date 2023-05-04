clear
clc
close all

% load the data

filename = 'subjData_mcmc_DeletingOutliers.mat';
load(filename);


predictedMPlus1_TW_BS_mean = nan(3,11,2);
predictedMPlus1_TW_BS_std = nan(3,11,2);
predictedMPlus1_MeanError_BS_mean = nan(3,11,2);
predictedMPlus1_MeanError_BS_std = nan(3,11,2);
predictedMPlus2_TW_BS_mean = nan(3,11);
predictedMPlus2_TW_BS_std = nan(3,11);
predictedMPlus2_MeanError_BS_mean = nan(3,11);
predictedMPlus2_MeanError_BS_std = nan(3,11);
 
for iExp = 1:11
    for iSub = 1:3
        if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_postcue)
            
            % multi-item trials
            
            % target weight for tne M-group experiements
            predictedMPlus1_TW_BS_mean(iSub,iExp,1) = mean(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF,'omitnan');
            predictedMPlus1_TW_BS_std(iSub,iExp,1) = std(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF,'omitnan');%/sqrt(length(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF));
            % mean error magnitude for tne M-group experiements
            predictedMPlus1_MeanError_BS_mean (iSub,iExp,1) =  mean(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError,'omitnan');
            predictedMPlus1_MeanError_BS_std (iSub,iExp,1) =  std(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError,'omitnan');%/sqrt(length(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError));
            
            % target weight for tne (M+1)-group experiements
            predictedMPlus1_TW_BS_mean(iSub,iExp,2) = mean(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF_plus1,'omitnan');
            predictedMPlus1_TW_BS_std(iSub,iExp,2) = std(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF_plus1,'omitnan')/sqrt(length(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF_plus1));
            
            % mean error magnitude for tne (M+1)-group experiements
            predictedMPlus1_MeanError_BS_mean(iSub,iExp,2) =  mean(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError_plus1,'omitnan');
            predictedMPlus1_MeanError_BS_std(iSub,iExp,2) =  std(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError_plus1,'omitnan')/sqrt(length(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError_plus1));
                       
            predictedMPlus2_TW_BS_mean(iSub,iExp) = mean(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF_plus2,'omitnan');
            predictedMPlus2_TW_BS_std(iSub,iExp) = std(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF_plus2,'omitnan')/sqrt(length(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanTF_plus2));
            

            predictedMPlus2_MeanError_BS_mean(iSub,iExp) =  mean(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError_plus2,'omitnan');
            predictedMPlus2_MeanError_BS_std(iSub,iExp) =  std(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError_plus2,'omitnan')/sqrt(length(subjData_mcmc(iSub).experiment(iExp).Bootstrap.meanError_plus2));

        end
    end
end
sprintf('done')

%%
toPlot1 = [2,3,7,8,9,10,11];
toPlot2 = [2,3,7,7,9,10];


toPlot3 = [2,3,3];
toPlot4 = [2,3,4];

color = [1,0.3250,0.0980;0.9290,0.6940,0.1250;0.3010,0.7450,0.9330];
marker_subj = {'o','o','o'};
marker_avg = {'-','-','-'};
interval = [-0.15,0,0.15];
h = figure();
h.Position = [-1910 -148 1382 1096];

pos = [];
for i = 1:6
 pos = [pos;[(mod(i-1,2))/2+(1/12) 1-(ceil(i/2))/3 5/12 1/3]];

end

markersize = 8;

stdCons = 2;
for isub = 1:3
    
    
    if isub~=1
        predictedMPlus1_TW_BS_mean(isub,toPlot2(end),2) = nan;
        predictedMPlus1_MeanError_BS_mean(isub,toPlot2(end),2) = nan;
    end

    %     subplot(3,2,isub*2-1)
    subplot('Position',pos(isub*2-1,:))

    h1 = errorbar([1:length(toPlot1)],predictedMPlus1_TW_BS_mean(isub,toPlot1,1),stdCons*predictedMPlus1_TW_BS_std(isub,toPlot1,1),marker_subj{isub},'color',color(3,:),'markersize',markersize,'LineWidth',1);hold on
    %     h2 = errorbar([2:length(toPlot1)]+interval(isub),predictedMPlus1_TW_BS_mean(isub,toPlot2,2),stdCons*predictedMPlus1_TW_BS_std(isub,toPlot2,2),marker_subj{isub},'color',color(1,:),'markersize',markersize,'LineWidth',1);hold on

    h1=  plot([1:length(toPlot1)],predictedMPlus1_TW_BS_mean(isub,toPlot1,1),'o-','color',color(3,:),'LineWidth',1,'markersize',markersize);hold on
    h2=  plot([2:length(toPlot1)],predictedMPlus1_TW_BS_mean(isub,toPlot2,2),'o-','color',color(1,:),'LineWidth',1,'markersize',markersize);hold on

    %     h3 = errorbar(toPlot3+1+interval(isub),predictedMPlus2_TW_BS_mean(isub,toPlot3), stdCons*predictedMPlus2_TW_BS_std(isub,toPlot3),marker_subj{isub},'color',color(2,:),'markersize',markersize,'LineWidth',1);hold on
    h3 = plot(toPlot4+1,predictedMPlus2_TW_BS_mean(isub,toPlot3),'o-','color',color(2,:),'LineWidth',1,'markersize',markersize);hold on
    %
    %   legend([h1,h2,h3],['S',num2str(isub),'''s data'],'M+1 prediction','M+2 prediction')

    xlim([0.5,length(toPlot1)+0.5])
    xticks([1:length(toPlot1)])
    ylabel('Target weight')
    ylim([0.3,1.1])
    %      xticklabels({'M3 (C3)';'  M4 \newline C3-->C4';'  M5 \newline C4-->C5 \newline C3-->C5';'  M6 \newline C5-->C6\newline C4-->C6';...
    %        '  M7\newline C3S3-->C4S3 ';'  M8\newline C4S3-->C4S4'})
%    xticklabels({'M3 (C3)';'  M4 \newline C3-->C4';'  M5 \newline C4-->C5 \newline C3-->C5';'  M6 \newline C5-->C6\newline C4-->C6';...
%        '  M7\newline C3S3-->C4S3 ';})
    set(gca,'fontSize',18)
%     set(gca,'XTick',[])
%     title(['S', num2str(isub), ': Multi-item trials'])    
   
%     subplot(3,2,isub*2)
     subplot('Position',pos(isub*2,:))
     h1 = errorbar([1:length(toPlot1)],predictedMPlus1_MeanError_BS_mean(isub,toPlot1,1),stdCons*predictedMPlus1_MeanError_BS_std(isub,toPlot1,1),marker_subj{isub},'color',color(3,:),'markersize',markersize,'LineWidth',1);hold on
%     h2= errorbar([2:length(toPlot1)]+interval(isub),predictedMPlus1_MeanError_BS_mean(isub,toPlot2,2),stdCons*predictedMPlus1_MeanError_BS_std(isub,toPlot2,2),marker_subj{isub},'color',color(1,:),'markersize',markersize,'LineWidth',1);hold on
    
    plot([1:length(toPlot1)],predictedMPlus1_MeanError_BS_mean(isub,toPlot1,1),'o-','color',color(3,:),'LineWidth',1,'markersize',markersize);hold on
    plot([2:length(toPlot1)],predictedMPlus1_MeanError_BS_mean(isub,toPlot2,2),'o-','color',color(1,:),'LineWidth',1,'markersize',markersize);hold on
    
%     h3 = errorbar(toPlot3+1+interval(isub),predictedMPlus2_MeanError_BS_mean(isub,toPlot3), stdCons*predictedMPlus2_MeanError_BS_std(isub,toPlot3),marker_subj{isub},'color',color(2,:),'markersize',markersize,'LineWidth',1);hold on
     plot(toPlot4+1,predictedMPlus2_MeanError_BS_mean(isub,toPlot3),'o-','color',color(2,:),'LineWidth',1,'markersize',markersize);hold on
%    
%     legend([h1,h2,h3],['S',num2str(isub),'''s data'],'M+1 prediction','M+2 prediction')
        
    xlim([0.5,length(toPlot1)+0.5])
    xticks([1:length(toPlot1)])
    ylabel('Mean Error Magnitude (pixels)')
%     xticklabels({'M3 (C3)';'  M4 \newline C3-->C4';'  M5 \newline C4-->C5 \newline C3-->C5';'  M6 \newline C5-->C6\newline C4-->C6';...
%        '  M7\newline C3S3-->C4S3 ';'  M8\newline C4S3-->C4S4'})
%     
    set(gca,'fontSize',18,'Ydir','reverse','FontName','Arial')
 ylim([10,120])
    
end
% %%
% position = [-1887 -399 984 1091];
