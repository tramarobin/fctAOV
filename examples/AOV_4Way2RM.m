clear
close all
clc

% functions path
getFunctionDirectory

% loading data and parameters
load projectParam
load vData

% savind directory
savedir=fullfile('4Way2RM');

% function
Tables=fctAOV(vData,projectParam.stats,projectParam.display,'Acceleration (m/s/s)',savedir);


