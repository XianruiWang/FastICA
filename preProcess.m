function [meanValue,whiteningMatrix,dewhiteningMatrix,processedSignal] = preProcess(inputSignal)
% Centering dataset and whitening it
% Usage: [meanValue,whiteningMatrix,dewhiteningMatrix,processedSignal] = preProcess(inputSignal)
% meanValue: mean of input signal
% whiteningMatrix: matrix transform original signal into uncorrelated signal
% dewhiteningMatrix: the inverse matrix of whiteningMatrix
% processedSignal: signal centered and whited
% inputSignal; original input signal to be processed, which is N*M matrix,
% N is the number of observations and M is the length of one observation 

% pcaNums; how many components will be reserved after PCA processing

% Centering: centeredSignal = inputSignal-mean(inputSignal)
% Whitening: processedSignal = ED^(-1/2)E'inputSignal where correlation
% PCA: discard unimportant component which can be seen as noise reduction and data compression
% matrix Rx=EDE'

% Some parameters
pcaNums = size(inputSignal,1);

% Centering
meanValue = mean(inputSignal,2);
centeredSignal = inputSignal-meanValue*ones(1,size(inputSignal,2));

% PCA
Rx = cov(centeredSignal', 1);
[E,D] = eig(Rx);
[D_valuesort,index] = sort(diag(D),'descend');
Rx_lenght=length(D_valuesort);
D_sort=zeros(Rx_lenght,Rx_lenght);
for i=1:Rx_lenght
    D_sort(i,i) = D_valuesort(i);
end 
E_sort = E(:,index);
pca_E(:,1:pcaNums) = E_sort(:,1:pcaNums);
pca_E(:,pcaNums+1:Rx_lenght)=E_sort(:,pcaNums+1:Rx_lenght);
pca_D(:,1:pcaNums) = D_sort(:,1:pcaNums);
pca_D(:,pcaNums+1:Rx_lenght)=D_sort(:,pcaNums+1:Rx_lenght);

% Whitening
whiteningMatrix = pca_E*sqrt(inv(pca_D))*pca_E';
dewhiteningMatrix = pca_E*sqrt(pca_D)*pca_E';
processedSignal = whiteningMatrix*centeredSignal;
fprintf('Data preprocessing accomplished \n');
