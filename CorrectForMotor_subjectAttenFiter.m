
figure(1)

folderName = '/Users/meow/Documents/MATLAB/PH.D in UCI/multiple centroids_exp data & analysis/Data Analysis_mcmc/';
filename = [folderName,'subjData_mcmc_DeletingOutliers.mat'];
load(filename);
%%

numGroups = [3,3,4,4,4,4,5,6,6,7,8];
numItemsPerGroups = [10,10,8,6,6,6,6,5,5,4,4];

items = numGroups.*numItemsPerGroups;

numSD = 2;

DiffPlus = cell(3,11,4);
DiffMinus = cell(3,11,4);

meanDiffPlus = nan(11,3,4);
meanDiffMinus = nan(11,3,4);

meanItems = nan(11,3,4);

for iSub = 1:3
    for iExp  = 1:11
        numColors = numGroups(iExp);
        N = numItemsPerGroups(iExp);
        totalItems = N*numColors;
        p = 0.01:0.01:0.99;
        dots = totalItems-totalItems*p;

        if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_precue)

            meanCurve_precue = subjData_mcmc(iSub).experiment(iExp).precue.subjEff.meanCurve;
            
            SD = sqrt(subjData_mcmc(iSub).experiment(iExp).precue.StandardError(:,1).^2+SEMotorError(iSub)^2);

            meanErrorPlusSTD = subjData_mcmc(iSub).experiment(iExp).precue.meanerror(:,1)-meanMotorError(iSub)+numSD*SD;
            meanErrorMinusSTD = subjData_mcmc(iSub).experiment(iExp).precue.meanerror(:,1)-meanMotorError(iSub)-numSD*SD;


            [~,idx_precue_plus] = min(abs(meanCurve_precue -  meanErrorPlusSTD'));
            numofTargetItemsPlus = dots(idx_precue_plus);

            [~,idx_precue_minus] = min(abs(meanCurve_precue - meanErrorMinusSTD'));
            numofTargetItemsMinus = dots(idx_precue_minus);

            [~,idx_precue] = min(abs(meanCurve_precue-(subjData_mcmc(iSub).experiment(iExp).precue.meanerror(:,1)-meanMotorError(iSub))'));
            numofTargetItems = dots(idx_precue);


            DiffPlus{iSub,iExp,2} = numofTargetItems - numofTargetItemsPlus;
            DiffMinus{iSub,iExp,2} = -numofTargetItems+numofTargetItemsMinus;

            meanDiffPlus(iExp,iSub,2) =  sqrt(sum(DiffPlus{iSub,iExp,2}.^2))/length(DiffPlus{iSub,iExp,2});
            meanDiffMinus(iExp,iSub,2) =  sqrt(sum(DiffMinus{iSub,iExp,2}.^2))/length(DiffMinus{iSub,iExp,2});

            meanItems(iExp,iSub,2) = mean(numofTargetItems);

        end


        if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_postcue)

            meanCurve_postcue = subjData_mcmc(iSub).experiment(iExp).postcue.subjEff.meanCurve;

            SD = sqrt(subjData_mcmc(iSub).experiment(iExp).postcue.StandardError(:,1).^2+SEMotorError(iSub)^2);

            meanErrorPlusSTD = subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,1)-meanMotorError(iSub)+numSD*SD;
            meanErrorMinusSTD = subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,1)-meanMotorError(iSub)-numSD*SD;


            [~,idx_postcue_plus] = min(abs(meanCurve_postcue -  meanErrorPlusSTD'));
            numofTargetItemsPlus = dots(idx_postcue_plus);

            [~,idx_postcue_minus] = min(abs(meanCurve_postcue - meanErrorMinusSTD'));
            numofTargetItemsMinus = dots(idx_postcue_minus);

            [~,idx_postcue] = min(abs(meanCurve_postcue-(subjData_mcmc(iSub).experiment(iExp).postcue.meanerror(:,1)-meanMotorError(iSub))'));
            numofTargetItems = dots(idx_postcue);


            DiffPlus{iSub,iExp,3} =  numofTargetItems- numofTargetItemsPlus;
            DiffMinus{iSub,iExp,3} = - numofTargetItems+numofTargetItemsMinus;

            meanDiffPlus(iExp,iSub,3) =  sqrt(sum(DiffPlus{iSub,iExp,3}.^2))/length(DiffPlus{iSub,iExp,3});
            meanDiffMinus(iExp,iSub,3) =  sqrt(sum(DiffMinus{iSub,iExp,3}.^2))/length(DiffMinus{iSub,iExp,3});

            meanItems(iExp,iSub,3) = mean(numofTargetItems);
        end

        if isfield(subjData_mcmc(iSub).experiment(iExp),'allCen') && ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_Onetarget)

            meanCurve_allcen = subjData_mcmc(iSub).experiment(iExp).allCen.subjEff.meanCurve;

            SD = sqrt(subjData_mcmc(iSub).experiment(iExp).allCen.StandardError(1)^2+SEMotorError(iSub)^2);

            meanErrorPlusSTD = subjData_mcmc(iSub).experiment(iExp).allCen.meanerror(1)-meanMotorError(iSub)+numSD*SD;
            meanErrorMinusSTD = subjData_mcmc(iSub).experiment(iExp).allCen.meanerror(1)-meanMotorError(iSub)-numSD*SD;




            [~,idx_allCen_plus] = min(abs(meanCurve_allcen-meanErrorPlusSTD));
            numofTargetItemsPlus = dots(idx_allCen_plus);

            [~,idx_allCen_minus] = min(abs(meanCurve_allcen-meanErrorMinusSTD));
            numofTargetItemsMinus = dots(idx_allCen_minus);


            [~,idx_allCen] = min(abs(meanCurve_allcen-(subjData_mcmc(iSub).experiment(iExp).allCen.meanerror(1)-meanMotorError(iSub))));
            numofTargetItems = dots(idx_allCen);


            DiffPlus{iSub,iExp,4} = numofTargetItems - numofTargetItemsPlus;
            DiffMinus{iSub,iExp,4} = -numofTargetItems+numofTargetItemsMinus;

            meanDiffPlus(iExp,iSub,4) =  sqrt(sum(DiffPlus{iSub,iExp,4}.^2))/length(DiffPlus{iSub,iExp,4});
            meanDiffMinus(iExp,iSub,4) =  sqrt(sum(DiffMinus{iSub,iExp,4}.^2))/length(DiffMinus{iSub,iExp,4});

            meanItems(iExp,iSub,4) = mean(numofTargetItems);
      
        end


        if ~isempty(subjData_mcmc(iSub).experiment(iExp).SummaryStats_Onetarget)
            
            meanCurve_Onetarget =  mean([subjData_mcmc(iSub).experiment(iExp).postcue.subjEff.meanCurve_perfect,subjData_mcmc(iSub).experiment(iExp).precue.subjEff.meanCurve_perfect],2);
            SD = sqrt(subjData_mcmc(iSub).experiment(iExp).Onetarget.StandardError(:,1).^2+SEMotorError(iSub)^2);

            meanErrorPlusSTD = subjData_mcmc(iSub).experiment(iExp).Onetarget.meanerror(:,1)-meanMotorError(iSub)+numSD*SD;
            meanErrorMinusSTD = subjData_mcmc(iSub).experiment(iExp).Onetarget.meanerror(:,1)-meanMotorError(iSub)-numSD*SD;

            [~,idx_noDis_plus] = min(abs(meanCurve_Onetarget- meanErrorPlusSTD'));
            numofTargetItemsPlus = dots(idx_noDis_plus);

            [~,idx_noDis_minus] = min(abs(meanCurve_Onetarget-meanErrorMinusSTD'));
            numofTargetItemsMinus = dots(idx_noDis_minus);


            [~,idx_noDis] = min(abs(meanCurve_Onetarget-(subjData_mcmc(iSub).experiment(iExp).Onetarget.meanerror(:,1)-meanMotorError(iSub))'));
            numofTargetItems = dots(idx_noDis);


            DiffPlus{iSub,iExp,1} = numofTargetItems - numofTargetItemsPlus;
            DiffMinus{iSub,iExp,1} = -numofTargetItems+numofTargetItemsMinus;

            meanDiffPlus(iExp,iSub,1) =  sqrt(sum(DiffPlus{iSub,iExp,1}.^2))/length(DiffPlus{iSub,iExp,1});
            meanDiffMinus(iExp,iSub,1) =  sqrt(sum(DiffMinus{iSub,iExp,1}.^2))/length(DiffMinus{iSub,iExp,1});

            meanItems(iExp,iSub,1) = mean(numofTargetItems);

        end
    end
end


order = [1:3,7:11,4:6];
upLimit = squeeze(sqrt(nansum(meanDiffPlus.^2,2))./sum(~isnan(meanDiffPlus),2));
lowerLimit = squeeze(sqrt(nansum(meanDiffMinus.^2,2))./sum(~isnan(meanDiffMinus),2));

upLimit = upLimit(order,:);
lowerLimit = lowerLimit(order,:);

meanItems_avg =squeeze(nanmean(meanItems(order,:,:),2));

interval = [-0.9:0.2:-0.1];

g1 = 1:8;
g2 = 9:11;

% color = [1,0.3250,0.0980;0.9290,0.6940,0.1250;0.3010,0.7450,0.9330;0.11,0.42,0.11];
color = [1,0,0;0.9290,0.6940,0.4250;0.3010,0.7450,0.9330;0,1,0];

totalItems = numGroups.*numItemsPerGroups;
tI = totalItems(order);

subplot(1,2,2)
edgealpha = 0.5;


for iType = 1:4
    p(iType) =plot(g1,meanItems_avg(g1,iType),'.-','color',color(iType,:),'LineWidth',1,'markersize',8);hold on
     
    if  iType==3
        y = [meanItems_avg(g1,iType)-upLimit(g1,iType); flipud(meanItems_avg(g1,iType)+lowerLimit(g1,iType))];
        h = fill([g1,fliplr(g1)]',y,color(iType,:));hold on
        h.EdgeAlpha = edgealpha;
        h.FaceAlpha = 0.05;
        h.EdgeColor = color(iType,:);
    end

    plot(g2,meanItems_avg(g2,iType),'.-','color',color(iType,:),'LineWidth',1,'markersize',8);hold on
    if  iType==3
        y = [meanItems_avg(g2,iType)-upLimit(g2,iType); flipud(meanItems_avg(g2,iType)+lowerLimit(g2,iType))];
        h = fill([g2,fliplr(g2)]',y,color(iType,:));hold on
        h.EdgeAlpha = edgealpha;
        h.FaceAlpha = 0.05;
        h.EdgeColor = color(iType,:);
    end

%     if iType~=1
        tt_g1  = mean(meanItems_avg(g1,iType),'omitnan');
        stdUP_g1 = mean(upLimit(g1,iType));
        stdLow_g1 = mean(lowerLimit(g1,iType));
        errorbar(g1(1)+interval(iType),tt_g1,stdUP_g1, stdLow_g1,'color',color(iType,:),'LineWidth',1,'marker','.','markersize',8);hold on

        tt_g2  = mean(meanItems_avg(g2,iType),'omitnan');
        stdUP_g2 = mean(upLimit(g2,iType));
        stdLow_g2 = mean(lowerLimit(g2,iType));
        errorbar(g2(end)-interval(iType),tt_g2,stdUP_g2, stdLow_g2,'color',color(iType,:),'LineWidth',1,'marker','.','markersize',8);hold on
%     end
end
plot(g1,tI(g1),'.-k'); hold on;
plot(g2,tI(g2),'.-k'); hold on;

xlim([-0.5,12.5])
ylim([0,34])

% legend([p(1),(2),p(3),p(4)],'zero-distractor','pre-cued','postcue','centroid-of-ALL')

xticks(1:11)
t = {'C3S0D10 NR \newline    3';'C3S0D10 \newline    3';'C4S0D8 \newline    4';'C4S0D6 \newline    4';'C2S2D6 \newline    4';'C0S4D6 \newline    4';...
    'C5S0D6 \newline    5' ;'C6S0D5 \newline    6 ';'C3S3D5 \newline    6';'C4S3D4 \newline    7';'C4S4D4 \newline    8'};
t = t(order);
xticklabels(t)
ylabel('number of stimulus items processed')
set(gca,'fontsize',16,'FontName','Arial')
set(0,'DefaultLineLineWidth',1)
title('Corrected for motor error subjects'' attention filters')