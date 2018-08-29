%% Palamedes toolbox
%% AG Mitchell 29.08.18
% Code running through the palamedes toolbox tutorial

%% Psychmetric functions(PFs)
% Choice of function
StimLevels = [1:1:6]; %data values for xaxis
pcorrect = PAL_Logistic([3 1 0.5 0], StimLevels);
% Generating crude plot of PF
figure(1)
plot(StimLevels, pcorrect, 'ko');

% Generating the same as above but with different parameter values
pcorrect = PAL_Logistic([3 1 0.5], StimLevels);
figure(2)
plot(StimLevels, pcorrect, 'ko');

% Maximum likeihood criterion for fitting PF, on a performance based 2AFC
% task 
StimLevels = [0.01 0.03 0.05 0.07 0.09 0.11];
NumPos = [45 55 72 85 91 100]; %number of trials in which observer gave correct response
OutOfNum = [100 100 100 100 100 100]; %number of trials for each stim level
% Specify type of function wish to use
PF = @PAL_Logistic; %logistic function
% Parameter values [threshold, slope, guess rate, lapse rate]
% Data around 75% correct when stim level around 0.05, high slope (more
% difficult to estimate), guess rate of 0.5 and assume lapse rate of 0
paramsValues = [0.05 50 0.5 0]; 
paramsFree = [1 1 0 0]; %two free parameters - thresh and slope, guess and lapse are fixed
% Curve fitting procedure
[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum,...
    paramsValues, paramsFree, PF); %PAL_PFML_Fit finds the best fitting parameters by way of interative 
%search through diff possible param values
% Graph showing data of smooth fitted function
PropCorrectData = NumPos./OutOfNum;
StimLevelsFine = [min(StimLevels):(max(StimLevels) - min(StimLevels))...
    ./1000:max(StimLevels)];
Fit = PF(paramsValues, StimLevelsFine);
% Plotting (proportion correct)
figure(3)
plot(StimLevels, PropCorrectData, 'k.', 'markersize', 40);
set(gca, 'fontsize', 12);
axis([0 .12 .4 1]);
hold on;
plot(StimLevelsFine, Fit, 'g-', 'linewidth', 4);

% Reached... estimating the errors