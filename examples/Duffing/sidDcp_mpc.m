% Script for testing the pqMPCm algorithm with the Duffing oscillator.
% pqMPCm : Model Predictive Control of a system based on a pqEDMD
% decompostion for Matlab.
% 
% The script assumes that the pqEDMDm package is also in the Matlab path.
%
% Load a pre-trained sidDecomposition object 
clear variables
load(fileparts(mfilename("fullpath"))+"/sidDcp.mat");

% 