% clear all
% This script plots S3 in the supplementary material
clc
folderName = '/Users/meow/Documents/MATLAB/PH.D in UCI/multiple centroids_exp data & analysis/Data Analysis_mcmc/';
filename = [folderName,'subjData_mcmc_DeletingOutliers.mat'];
load(filename);
%%
numGroups = [3,3,4,4,4,4,5,6,6,7,8];
numItemsPerGroups = [10,10,8,6,6,6,6,5,5,4,4];

p = 0.01:0.01:0.99;

expIdx = [2,3,7:10];
colors = rand(3,7);

meanMotorError =[14.7423,14.0099, 13.0468];
markersize = 10;

figure()
for i = 1:length(expIdx)

   
   iExp = expIdx (i);

    numColors = numGroups(iExp);
    N = numItemsPerGroups(iExp);
    totalItems = N*numColors;

    for iSub = 1:3
        subjectNum = iSub;
        
        if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_precue) 
            
            subplot(3,6,(iSub-1)*6+i)
            dots = totalItems-totalItems*p;


            meanCurve = mean(subjData_mcmc(iSub).experiment(iExp).postcue.subjEff.meanCurve_perfect,2);

            meanError = mean(subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,1));
            
            [~,idx] = min(abs(meanCurve-meanError));
            numofTargetItems = dots(idx);
             

            [~,idx_motor] = min(abs(meanCurve-(meanError-meanMotorError(iSub))));
            numofTargetItems_motor = dots(idx_motor);


            meanError_precue = mean(subjData_mcmc(iSub).experiment(iExp).precue.meanerror(:,1));
            
            [~,idx] = min(abs(meanCurve-meanError_precue));
            numofTargetItems_precue = dots(idx);
             

            [~,idx_motor] = min(abs(meanCurve-(meanError_precue-meanMotorError(iSub))));
            numofTargetItems_motor_precue = dots(idx_motor);
         

            plot(meanCurve,dots,'k--'); hold on

            h1 = plot(meanError,numofTargetItems,'r*','MarkerSize',markersize); hold on
            h2 = plot(meanError-meanMotorError(iSub),numofTargetItems_motor,'b*','MarkerSize',markersize); hold on           
            h3 = plot(meanError_precue,numofTargetItems_precue,'r+','MarkerSize',markersize); hold on
            h4 = plot(meanError_precue-meanMotorError(iSub),numofTargetItems_motor_precue,'b+','MarkerSize',markersize); hold on
           
 
            ylim([ 0 32])
            xlim([0 180])
            title(['M = ', num2str(numColors)])
%             xlabel('mean error mag')
%             ylabel('number of items processed')
            set(gca,'fontsize',14)
            
        
        end

    end

end
  legend([h1;h2;h3;h4],{'post-cued:uncorrect for motor'; ' post-cued:corrected for motor';'pre-cued:uncorrect for motor';' pre-cued:corrected for motor'})
            

%%

