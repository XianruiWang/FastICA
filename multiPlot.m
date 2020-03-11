function multiPlot(signalName,signal)
figure('NumberTitle', 'off', 'Name', signalName);
figureNums = size(signal,1);
set(gcf,'unit','centimeters','position',[20 5 25 12]);
for figureIter = 1:figureNums
    subplot(figureNums,1,figureIter);
    plot(signal(figureIter,:));
end