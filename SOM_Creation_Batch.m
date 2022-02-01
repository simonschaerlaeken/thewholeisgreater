% ---------------------------- %
%          SOM Batch
% ---------------------------- %
clear *;
clc;
close all;
%% Read data
table_raw = readtable('/YourPath/Table_Data_ZScored_Residuals.csv');
[num, numLabel, text, textLabel] = convert_data_NA(table_raw);
% Change _ in .
numLabel = strrep(numLabel,'_', '.');
%% Select Data
selectVariables = ...
{'zscore.TNT.dprime.tap.trace.residuals','zscore.SAAT.rt.sd.correct.sustained.residuals','zscore.SAAT.dprime.impulsive.residuals','zscore.FLANKER.rcs.overall.residuals','zscore.BOXED.rt.mean.correct.residuals','zscore.TASKSWITCH.rcs.overall.residuals','zscore.STROOP.rcs.overall.residuals','zscore.SPATIALSPAN.object.count.span.residuals','zscore.BACKWARDSSPATIALSPAN.object.count.span.residuals'};

% Find index selected variables
indexSelectVariables = zeros(numel(selectVariables),1);
for i = 1:numel(selectVariables)
    variableName = selectVariables{i};
    for j = 1:numel(numLabel)
        columnName = numLabel{j};
        if strcmp(variableName,columnName)
            indexSelectVariables(i) = j;
            break
        end
    end
end
data = num(:,indexSelectVariables);
dataLabel = numLabel(:,indexSelectVariables);

%% Create structure & Normalize
sD = som_data_struct(data,'name','Test data',...
                                 'comp_names',dataLabel);

%% Batch training
% Configs
initMethods = {'randinit', 'lininit'};
trainMethods = {'seq','batch'};
neighMethods = {'gaussian','cutgauss','bubble','ep'};
nbIter = 10;
nbIterTest = 10;
thresholdTest = 0.7;

% Lattice Size -- Kohonen: 5?n nodes where n is the number of data point
% Kohonen also suggests that the optimal ratio of height to width of
% the lattice is equal to the ratio of the two largest eigenvalues of the 
% autocorrelation matrix
% sMTable = struct('type',{},...
%                  'codebook',{},...
%                  'topol',{},...
%                  'labels',{},...
%                  'neigh',{},...
%                  'mask',{},...
%                  'trainhist',{},...
%                  'name',{},...
%                  'comp_names',{},...
%                  'comp_norm',{});
errorTable = []; 
lowestQeTe = 9999;
lowestQeTeSM = [];
lowestTe = 9999;
lowestTeSM = [];
lowestQe = 9999;
lowestQeSM = [];
lowestCbe = 9999;
lowestCbeSM = [];

% Batch training
for iInitMethods = 1:numel(initMethods)
    for iTrainMethods = 1:numel(trainMethods)
        for iNeighMethods = 1:numel(neighMethods)
            for iNbIter = 1:nbIter
                for iNbIterTest = 1:nbIterTest
                    % Print
                    disp(['---------',initMethods{iInitMethods},'---',...
                        trainMethods{iTrainMethods},'---',neighMethods{iNeighMethods},'---',...
                        num2str(iNbIter), '--- test: ',num2str(iNbIterTest)])
                    
                    % Divide the data into training set and testing set
                    [m,~] = size(sD.data);
                    idx = randperm(m);
                    trainData = sD.data(idx(1:round(thresholdTest*m)),:); 
                    testData = sD.data(idx(round(thresholdTest*m)+1:end),:) ;

                    % Create structure & Normalize
                    sDtrain = som_data_struct(trainData,'name','Test data',...
                                 'comp_names',dataLabel);
                    sDtest = som_data_struct(testData,'name','Test data',...
                                 'comp_names',dataLabel);

                    % Create map
                    sM = som_make(sDtrain,...
                        'init', initMethods{iInitMethods},...
                        'algorithm', trainMethods{iTrainMethods},...
                        'neigh', neighMethods{iNeighMethods});
                    [qeTrain,teTrain,cbeTrain] = som_quality(sM,sDtrain);
                    [qeTest,teTest,cbeTest] = som_quality(sM,sDtest);      
                    
                    % Test for QE & TE
                    errorQeTe = mean([real(qeTrain),real(teTrain),real(qeTest),real(teTest)]);
                    if errorQeTe<lowestQeTe
                        lowestQeTe = errorQeTe;
                        lowestQeTeSM = sM;
                    end
                    % Test for TE
                    errorTe = mean([real(teTrain),real(teTest)]);
                    if errorTe<lowestTe
                        lowestTe = errorTe;
                        lowestTeSM = sM;
                    end
                    % Test for QE
                    errorQe = mean([real(qeTrain),real(qeTest)]);
                    if errorQe<lowestQe
                        lowestQe = errorQe;
                        lowestQeSM = sM;
                    end
                    % Test for CBE
                    errorCbe = mean([real(cbeTrain),real(cbeTest)]);
                    if errorCbe<lowestCbe
                        lowestCbe = errorCbe;
                        lowestCbeSM = sM;
                    end
                    errorTable(end+1,:) = [errorQeTe,errorTe,errorQe,errorCbe,real(qeTrain),real(teTrain),real(cbeTrain),real(qeTest),real(teTest),real(cbeTest)];
                end
            end
        end
    end    
end

% Save
save('iLead_sMTable.mat','sMTable')
save('iLead_errorTable.mat','errorTable')

% Find index
indexLowestQeTe = find(errorTable(:,1)==lowestQeTe);
indexLowestQe = find(errorTable(:,2)==lowestQe);
indexLowestCbe = find(errorTable(:,3)==lowestCbe);

% Choose best model
clear sM;
disp(['Minimum Error QETE: ', num2str(lowestQeTe), ' - ', num2str(indexLowestQeTe)]);
disp(['Minimum Error QE: ', num2str(lowestQe), ' - ', num2str(indexLowestQe)]);
disp(['Minimum Error CBE: ', num2str(lowestCbe), ' - ', num2str(indexLowestCbe)]);
colorTable = linspecer(size(errorTable,2));
% Visualise with figure
figure;
hold on;
for i = 1:size(errorTable,2)
    plot(errorTable(:,i),'lineWidth',2, 'Color', colorTable(i,:));
end
vline(indexLowestQeTe,'g',['Minimum Error QETE: ', num2str(lowestQeTe)])
vline(indexLowestQe,'b',['Minimum Error QE: ', num2str(lowestQe)])
vline(indexLowestCbe,'r',['Minimum Error CBE: ', num2str(lowestCbe)])
hold off;
legend({'errorQeTe','errorQe','errorCbe','qeTrain','teTrain','cbeTrain','qeTest','teTest','cbeTest'})
set(gcf,'Name','ErrorTable SOM')

% Save best models
sM = lowestQeTeSM;
save('iLead_sM_QETE.mat','sM')
sM = lowestQeSM;
save('iLead_sM_QE.mat','sM')
sM = lowestCbeSM;
save('iLead_sM_CBE.mat','sM')

load('iLead_sM_CBE.mat')
%% TESTING ON DIFFERENT COHORT and TIME
% Verify that our model doesn't differ based on cohort and time

% Select in Dataset
dataCohort = num(:,1);
idx_C3 = dataCohort==3;
idx_C5 = dataCohort==5;
idx_C7 = dataCohort==7;
dataTimeCell = text(:,2);
dataTime = zeros(size(dataTimeCell));
for i = 1:size(dataTimeCell,1)
    dataTime(i) = str2num(dataTimeCell{i}(2));
end
idx_T1 = dataTime==1;
idx_T2 = dataTime==2;
idx_T3 = dataTime==3;
idx_T4 = dataTime==4;

% Create structures
testData_C3 = data(idx_C3,:);
sDtest_C3 = som_data_struct(testData_C3,'name','Test data',...
                                 'comp_names',dataLabel);
testData_C5 = data(idx_C5,:);
sDtest_C5 = som_data_struct(testData_C5,'name','Test data',...
                                 'comp_names',dataLabel);
testData_C7 = data(idx_C7,:);
sDtest_C7 = som_data_struct(testData_C7,'name','Test data',...
                                 'comp_names',dataLabel);
testData_T1 = data(idx_T1,:);
sDtest_T1 = som_data_struct(testData_T1,'name','Test data',...
                                 'comp_names',dataLabel);
testData_T2 = data(idx_T2,:);
sDtest_T2 = som_data_struct(testData_T2,'name','Test data',...
                                 'comp_names',dataLabel);
testData_T3 = data(idx_T3,:);
sDtest_T3 = som_data_struct(testData_T3,'name','Test data',...
                                 'comp_names',dataLabel);
testData_T4 = data(idx_T4,:);
sDtest_T4 = som_data_struct(testData_T4,'name','Test data',...
                                 'comp_names',dataLabel);
% Test Error
cbeTest = zeros(7,1);
[~,~,cbeTest(1)] = som_quality(sM,sDtest_C3);
[~,~,cbeTest(2)] = som_quality(sM,sDtest_C5);
[~,~,cbeTest(3)] = som_quality(sM,sDtest_C7);
[~,~,cbeTest(4)] = som_quality(sM,sDtest_T1);
[~,~,cbeTest(5)] = som_quality(sM,sDtest_T2);
[~,~,cbeTest(6)] = som_quality(sM,sDtest_T3);
[~,~,cbeTest(7)] = som_quality(sM,sDtest_T4);
% Save
save('iLead_cbe_Cohort_Time.mat','cbeTest')

%% Color map
% Create the yellow-grey color map
n = 20;    %// number of colors

R1 = linspace(4/255,65/255,n);  %// Red from 1 to 0
G1 = linspace(4/255,65/255,n);   %// Green all zero
B1 = linspace(4/255,65/255,n);  %// Blue from 0 to 1

R2 = linspace(65/255,251/255,n);  %// Red from 1 to 0
G2 = linspace(65/255,188/255,n);   %// Green all zero
B2 = linspace(65/255,4/255,n);  %// Blue from 0 to 1

R3 = linspace(251/255,255/255,n);  %// Red from 1 to 0
G3 = linspace(188/255,235/255,n);   %// Green all zero
B3 = linspace(4/255,185/255,n);  %// Blue from 0 to 1

R = [R1,R2,R3];
G = [G1,G2,G3];
B = [B1,B2,B3];
cmap = colormap( [R(:), G(:), B(:)] );  %// create colormap


%% Plot Individual feature maps
figure;som_show(sM);
for iComp = 1:numel(sM.comp_names)
    figure;som_show(sM, 'comp', iComp,'colormap',cmap);
    set(gcf,'NumberTitle','off') %don't show the figure number
    set(gcf,'Name',sM.comp_names{iComp})
end

% Plot hitmap
sHits = som_hits(sM,sD);
figure;som_show(sM,'empty','Labels');
som_show_add('hit',sHits);
set(gcf,'Name','Hits')
%% k-mean Clustering
[c, p, err, ind] = kmeans_clusters(sM);

% Error + Size Cluster and Silhouette
[minError, minIndex] = min(err);
figure; plot(err,'lineWidth',2);title('Error KMeans');
set(gcf,'Name','ScreePlot Error sMap')
figure; plot(ind,'lineWidth',2);title('Davies-Bouldin index');
set(gcf,'Name','ScreePlot Davies-Bouldin index sMap')
 
% Find number of participants for each clusters and silhouette
bmu = som_bmus(sM,sD);
indexP = 2:length(p);
silhouetteTable = zeros(size(indexP));
sizeClusters = zeros(length(p));
for i = 1:length(p)
    particClusters = zeros(size(bmu,1),1);
    for iPart = 1:size(bmu,1)
        particClusters(iPart) = p{i}(bmu(iPart));
    end
    uniqueClusters = unique(particClusters);
    for j = 1:length(uniqueClusters);
        sizeClusters(i,j) = sum(particClusters==uniqueClusters(j));
    end
    % Silhouette
    if i ~= 1
        silhouetteTable(i-1) = mean(silhouette(sM.codebook,p{i},'Euclidean'));
    end
end

% Plot Silhouette table
figure; plot(indexP,silhouetteTable, 'lineWidth',2);title('Silhouette');
set(gcf,'Name','Mean Silhouette')

% Disp and Save size Cluster
disp(sizeClusters)
T_sizeCluster = table(sizeClusters);
writetable(T_sizeCluster,'iLead_sizeCluster.csv');

%% Choose nbCluster
nCluster = 6; % Choose based on graphs printed before
pSel = p{nCluster};
figure;silhouette(sM.codebook,p{nCluster},'Euclidean');
set(gcf,'Name','Silhouette')

% If ReOrder -- to rank profiles from worse to best
if 1
    oldOrder = [1 2 3 4 5 6];
    newOrder = [1 2 5 4 6 3];
    pSelOrdered = change_order_cluster( pSel,  oldOrder, newOrder );
    pSel = pSelOrdered;
end

% Save new profile assignment
save('iLead_pSel.mat','pSel')

% Create Color map
colorList = [
    102,194,165; 
    252,141,98;         
    141,160,203;
    231,138,195;
    166,216,84;
    254,216,68; 
    229,196,148;
    179,179,179]/255; 

colorTable = zeros(length(pSel),3);
for i = 1:length(pSel)
    colorTable(i,:) = colorList(pSel(i),:);
end

% Print profile map
figure;som_cplane(sM,colorTable);
set(gcf,'Name','Color KMean Map')

%% Find association cluster
bmu = som_bmus(sM,sD);
particClusters = [];
for iPart = 1:size(bmu,1)
    particClusters(iPart,1) = pSel(bmu(iPart));
    particClusters(iPart,2) = bmu(iPart);
end

% Write table
T_particClusters = table(num,particClusters);
writetable(T_particClusters,'iLead_ParticCluster.csv');
save('iLead_ParticCluster.mat','particClusters')
%load('iLead_ParticCluster.mat')
T_Clusters = table(sM.codebook,pSel);
writetable(T_Clusters,'iLead_Clusters.csv');

%% Quick comparisons between EFs (T-TEST)
dataStat = {};
dataStatAvg = [];
dataStatStd = [];
dataLegend = {};

sDraw = sD;
%sDraw = som_denormalize(sD);
for iClust = 1:nCluster
    dataClust = sDraw.data(particClusters(:,1)==iClust,:);
    dataStat{iClust} = dataClust;
    dataStatAvg(iClust,:) = nanmean(dataClust,1);
    dataStatStd(iClust,:) = nanstd(dataClust,1);
    dataLegend{iClust} = strcat('Cluster ', num2str(iClust)); 
end

% Draw T Test X - LINE
% Figure
figure
hold on
colorMatrix = colorList;
xAbs = 1:size(dataStatAvg,2);
% For each line in the matrix
for i = 1:size(dataStatAvg,1)
    meanValues = dataStatAvg(i,:); % tes valeurs moyennes
    highIC = dataStatStd(i,:); %les valeurs d'intervalle de confiance haut 
    lowIC = dataStatStd(i,:);%les valeurs d'intervalle de confiance bas 
    colorVector = colorMatrix(i,:);
    errorbar(xAbs,meanValues,lowIC,highIC,'Color', colorVector ,'MarkerSize',25,'linewidth',1.5);
end
legend(dataLegend)
xTicks = sD.comp_names;               
set(gca, 'XTick', 1:numel(xTicks), 'XTickLabel', xTicks)
set(gca,'XTickLabelRotation',45)
set(gcf,'Name','DistributionCluster')

% Draw T Test X - BAR
dataStatAvgCondensed = zeros(size(dataStatAvg,1),3);
dataStatAvgCondensed(:,1) = mean(dataStatAvg(:,1:3),2); % AC
dataStatAvgCondensed(:,2) = mean(dataStatAvg(:,4:7),2); % IR
dataStatAvgCondensed(:,3) = mean(dataStatAvg(:,8:9),2); % WM
%dataStatAvgCondensed(:,3) = mean(dataStatAvg(:,8:11),2); % WM

dataStatStdCondensed = zeros(size(dataStatStd,1),3);
dataStatStdCondensed(:,1) = mean(dataStatStd(:,1:3),2); % AC
dataStatStdCondensed(:,2) = mean(dataStatStd(:,4:7),2); % IR
dataStatStdCondensed(:,3) = mean(dataStatStd(:,8:9),2); % WM
%dataStatStdCondensed(:,3) = mean(dataStatStd(:,8:11),2); % WM

% Plot Figure
label = {'AC','IR', 'WM'};
tickLabel = {};
figure()
hold on
idx = 1;
for i = 1:size(dataStatAvgCondensed,1)
    for j = 1:size(dataStatAvgCondensed,2)
        meanValues = dataStatAvgCondensed(i,j); % tes valeurs moyennes
        highIC = dataStatStdCondensed(i,j); %les valeurs d'intervalle de confiance haut 
        lowIC = dataStatStdCondensed(i,j);%les valeurs d'intervalle de confiance bas 
        colorVector = colorMatrix(i,:);
        h=bar(idx,meanValues);
        set(h,'FaceColor',colorVector);
        errorbar(idx,meanValues,lowIC,highIC,'Color', colorVector ,'MarkerSize',25,'linewidth',1.5);
        tickLabel{end+1} = [num2str(i),'.',label{j}];
        idx = idx + 1;
    end
end
hold off
set(gca, 'XTick', 1:idx-1, 'XTickLabel', tickLabel)
set(gcf,'Name','DistributionClusterBar')


%% Quick comparisons between profiles for Math and Reading (T-TEST)
selectVariables = ...
{'zscore.MATH.FLUENCY.acc.sum.overall','zscore.READING.COMPREHENSION.acc.score.overall'};

indexSelectVariables = zeros(numel(selectVariables),1);
for i = 1:numel(selectVariables)
    variableName = selectVariables{i};
    for j = 1:numel(numLabel)
        columnName = numLabel{j};
        if strcmp(variableName,columnName)
            indexSelectVariables(i) = j;
            break
        end
    end
end
dataY = num(:,indexSelectVariables);
dataYLabel = numLabel(:,indexSelectVariables);

% T-TEST Y
dataYStat = {};
dataYStatAvg = [];
dataYStatStd = [];
dataYLegend = {};

for iClust = 1:nCluster
    dataClust = dataY(particClusters(:,1)==iClust,:);
    dataYStat{iClust} = dataClust;
    dataYStatAvg(iClust,:) = nanmean(dataClust,1);
    dataYStatStd(iClust,:) = nanstd(dataClust,1);
    dataYLegend{iClust} = strcat('Cluster ', num2str(iClust)); 
end

% PLOT - BAR graph
label = {'Math','Reading'};
tickLabel = {};
figure()
hold on
idx = 1;
for i = 1:size(dataYStatAvg,1)
    for j = 1:size(dataYStatAvg,2)
        meanValues = dataYStatAvg(i,j); % tes valeurs moyennes
        highIC = dataYStatStd(i,j); %les valeurs d'intervalle de confiance haut 
        lowIC = dataYStatStd(i,j);%les valeurs d'intervalle de confiance bas 
        colorVector = colorMatrix(i,:);
        h=bar(idx,meanValues);
        set(h,'FaceColor',colorVector);
        errorbar(idx,meanValues,lowIC,highIC,'Color', colorVector ,'MarkerSize',25,'linewidth',1.5);
        tickLabel{end+1} = [num2str(i),'.',label{j}];
        idx = idx + 1;
    end
end
hold off
set(gca, 'XTick', 1:idx-1, 'XTickLabel', tickLabel)
set(gcf,'Name','MathReading Cluster')
% ---------------------------------------------------------