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
[paramsValues, LL, exitflag] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum,...
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

% Estimating the errors
% Bootstrapping function - simulates actual experiment repeatedly to
% estimate SD of threshold and slope. Requires PF fitting routine to have
% already been run
B = 400; %how much simulated data estimated
[SD, paramsSim, LLSim, converged] = PAL_PFML_BootstrapParametric...
    (StimLevels, OutOfNum, paramsValues, paramsFree, B, PF); %bootstrap function, need to use the same paramsValues as estimated with PAL_PFML_Fit
% SD output is the estimate of the errors for: [thresh slope guess lapse]

% Estimating the goodness of fit
% Determined by comparing the two models statistically
% Outputs: Dev = deviance and pDev = goodness of fit measure (always
% between 0 and 1), uses the same args as the error estimation routine
B = 1000;
[Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, ...
    OutOfNum, paramsValues, paramsFree, B, PF);

%% Putting it all together

