%%% The main script that calls other functions to reproduce the figures in:
%%% Systems biology reveals uncoupling beyond UCP1 in human white fat-derived beige adipocytes
%%% Author: Elin Nyman, Linköping University

%% clear workspace before running the script
close all
clear variables

%% globals are used by all functions and scripts
global BartData
global BartData2
global BartData4
global pNamesOpt  
global icOrig 
global modelName
global maxsim
global minsim
global maxsim2
global minsim2
global maxsim3
global minsim3
global maxsim4
global minsim4
global maxsim5
global minsim5

%% load data from the folder DATA
dd=cd;
cd DATA
load('Bart.mat')
load('BartData2.mat')
load('BartData4.mat')
cd(dd)
BartData = Bart;

%% load collected acceptable parameter values
load('parameterValues.dat')

%% mex compile the model file and define model parameters
modelName = 'model_H3'; %Replace with H1, H2 or H3 to change model
optModel = IQMmodel(strcat(modelName,'.txt'));
tmpModel = optModel;
IQMmakeMEXmodel(optModel,modelName);
[pNamesOpt, pValues] = IQMparameters(optModel);
icOrig = IQMinitialconditions(optModel);

%% decide parameters to plot based on number of data points, i.e. degrees of freedom in chi2 test
limit=chi2inv(0.95,10*4+9*5);
index=find(parameterValues(:,1)<limit);
goodValues=parameterValues(index,:);
sizeparam = size(goodValues(1,:),2);
maxparam = zeros(sizeparam,1);
minparam = zeros(sizeparam,1);
indexmin = zeros(sizeparam,1);
indexmax = zeros(sizeparam,1);
minparams = zeros(sizeparam);
maxparams = zeros(sizeparam);

for i = 1:sizeparam
    minparam(i)=min(goodValues(:,i));
    maxparam(i)=max(goodValues(:,i),[],1);
    indexmin(i)=find(goodValues(:,i)<=minparam(i),1);
    indexmax(i)=find(goodValues(:,i)>=maxparam(i),1);
    minparams(i,:)=goodValues(indexmin(i),:);
    maxparams(i,:)=goodValues(indexmax(i),:);
end
uniqueparams=unique(goodValues,'rows');
sizeunique=size(uniqueparams,1);

maxsim=zeros(6,73);
minsim=1000.*ones(6,73);
maxsim2=zeros(6,55);
minsim2=1000.*ones(6,55);
maxsim3=zeros(6,55);
minsim3=1000.*ones(6,55);
maxsim4=zeros(6,73);
minsim4=1000.*ones(6,73);
maxsim5=zeros(6,73);
minsim5=1000.*ones(6,73);

%% call costfunction to simulate the model for the given parameters and collect the range of possible model behaviors
for i = 2:100:sizeunique
     costfunction(uniqueparams(i,:),0);
end

%% plot the reults from the full range of possible model behaviors as gray areas
plotscript(minsim,minsim2,minsim3,minsim4,minsim5,maxsim,maxsim2,maxsim3,maxsim4,maxsim5)

%% plot the best solution as lines
costfunction(minparams(1,:),1);