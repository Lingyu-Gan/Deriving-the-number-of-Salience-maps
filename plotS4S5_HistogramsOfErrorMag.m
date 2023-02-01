% This script is to plot S4 and S5 in the supplementary material

clc
close all
clear 
folderName = '/Users/meow/Documents/MATLAB/PH.D in UCI/multiple centroids_exp data & analysis/Data Analysis_mcmc/';
filename = [folderName,'subjData_mcmc_DeletingOutliers.mat'];
load(filename);

numGroups = [3,3,4,4,4,4,5,6,6,7,8];
numItemsPerGroups = [10,10,8,6,6,6,6,5,5,4,4];

expIdx = [2,3,7:10];
colors = rand(3,7);

for i = 1:length(expIdx)

    iExp = expIdx (i);

    numColors = numGroups(iExp);
    N = numItemsPerGroups(iExp);
    totalItems = N*numColors;

    edges = [0:20:320];% x-axis

    for iSub = 1:3
        subjectNum = iSub;
         
         errorVec = [];
         errorVec_s = [];

        if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_precue)

            for iColor = 1:numColors 
           

                answer_s = [subjData_mcmc(iSub).experiment(iExp).postcue.ItemPosAndResp_s{iColor,1}(:,iColor),...
                    subjData_mcmc(iSub).experiment(iExp).postcue.ItemPosAndResp_s{iColor,2}(:,iColor)];

                resp_s = subjData_mcmc(iSub).experiment(iExp).postcue.ItemPosAndResp_s{iColor,3};

                errorVec_s = [errorVec_s;sqrt(sum((resp_s-answer_s ).^2,2))];


                answer = [subjData_mcmc(iSub).experiment(iExp).postcue.ItemPosAndResp{iColor,1}(:,iColor),...
                    subjData_mcmc(iSub).experiment(iExp).postcue.ItemPosAndResp{iColor,2}(:,iColor)];

                resp = subjData_mcmc(iSub).experiment(iExp).postcue.ItemPosAndResp{iColor,3};

                errorVec = [errorVec;sqrt(sum((resp-answer).^2,2))];


            end
           
           figure(1) % for multi-item trials
           subplot(3,6,i+(iSub-1)*6)
           histogram(errorVec,edges,'Normalization','probability');
           text(200,0.4, num2str(round(mean(errorVec),2)))
           ylim([0,0.6])
           title(['multi M = ', num2str(numColors)])
%            xlabel('error magnitude')
%            ylabel('prob')
           set(gca,'fontsize',14)
           
           figure(2) % for singletons trials
           subplot(3,6,i+(iSub-1)*6)
           histogram(errorVec_s,edges,'Normalization','probability');
           text(200,0.4, num2str(round(mean(errorVec_s),2)))
           ylim([0,0.6])
           title(['singleton M = ', num2str(numColors)])
%            xlabel('error magnitude')
%            ylabel('prob')
           set(gca,'fontsize',14)




        end

    end

end
