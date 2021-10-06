% Fast ICA algorithm
% A: the mixing matrix
% W: the hypothesised unmixing matrix
% B: the space spanned by unmixing vector
% w: unmixing vector 
% wOld: w in last iteration
% failNums: fow many times will try until give up
% round: the number of Independent Component under searching
% icNums: total numbers of Independent Components under searching
% iteration: iteration numbers for searching one Independent Component
% Author: Xianrui Wang, from Center of Intellegent Acoustics and Immersive 
% Communications(CIAIC), copyright reserved
% Contact me: wangxianrui@mail.nwpu.edu.cn
% For more resource of Blind Source Separation, visit wangxianrui.club

%%
%load data and precess it
load('origSignal');
% preprocessing
inputSignal=randn(size(origSignal,1))*origSignal;
[meanValue,whiteningMatrix,dewhiteningMatrix,processedSignal] = preProcess(inputSignal);
% some parameters
icNums = size(processedSignal,1);
numSamples = size(processedSignal,2);
A = zeros(icNums);
B = zeros(icNums);
W = zeros(icNums);
failNums = 2000;
epsilon=0.0001;
a1=1;
round=1;
%%
% fixed point algorithm for fast ICA 
while round <= icNums
    w = randn(icNums,1);
    wOld = zeros(icNums,1);
    iteration = 1;
    while iteration<=failNums
        iteration=iteration+1;
        w = w-B*B'*w;
        w = w/norm(w);
        if(iteration==failNums)
            fprintf('algorithm not converge this time \n');
        end
        if(norm(w - wOld) < epsilon || norm(w + wOld) < epsilon)
            fprintf('algorithm has converged \n');
            B(:,round) = w;
            A(:,round) = dewhiteningMatrix * w;
            W(round,:) = w' * whiteningMatrix;
            round= round+1;
            break;
        else
            wOld = w;
            hypTan = myTanh(a1 * processedSignal' * w);
            w = (processedSignal * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / numSamples;
        end
    end
end

outputSignal = W*inputSignal+A\meanValue;
%%
% figures
multiPlot('original signal',origSignal);
multiPlot('mixed signal',inputSignal);
multiPlot('output siganl',outputSignal);


