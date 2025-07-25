clear all
close all
clc

% base RDM matrix (2 by 2)
baseRDM = [0 1; 1 0];

% reward coding matrix
X1 = squareform(round(imresize(baseRDM,8)))';

% task coding matrix
X2 = squareform(repmat(round(imresize(baseRDM,4)),2,2))';

% relevant/irrelevant features coding matrix
corner1 = repmat(baseRDM,2,2);
corner2 = round(imresize(baseRDM,2));
X3 = squareform(repmat([corner1 ones(4); ones(4) corner2],2,2))'; % relevant
X4 = squareform(repmat([corner2 ones(4); ones(4) corner1],2,2))'; % irrelevant

% motor coding matrix
corner3 = round(imresize(repmat(baseRDM,1,2),[4 4]));
X5 = squareform(repmat([corner1 corner3'; corner3 corner2],2,2))';

X = zscore([X1 X2 X3 X4 X5]);

%%

model1 = fitlm(X(:,[2 3 4 5]),X(:,1));
model2 = fitlm(X(:,[1 3 4 5]),X(:,2));
model3 = fitlm(X(:,[1 2 4 5]),X(:,3));
model4 = fitlm(X(:,[1 2 3 5]),X(:,4));
model5 = fitlm(X(:,[1 2 3 4]),X(:,5));

VIF1 = 1/(1-model1.Rsquared.Ordinary)
VIF2 = 1/(1-model2.Rsquared.Ordinary)
VIF3 = 1/(1-model3.Rsquared.Ordinary)
VIF4 = 1/(1-model4.Rsquared.Ordinary)
VIF5 = 1/(1-model5.Rsquared.Ordinary)