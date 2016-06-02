function [SEIR]= mrQ_fitSEIR_T1(SEIRdir,outDir,checkData)
% 
%  [SEIR]=mrQ_fitSEIR_T1_ver2(SEIRdir,outDir,[checkData=1])
% 
% Loads SEIR_Dat.mat file from the directory 'data' in the SEIRdir, and
% performs T1 mapping. 
% 
% 
% INPUTS:
%   SEIRdir   -   Path to your desired T1 image. Use the same path as in
%                 getSEIR.
% 
%   checkData -   If you want to visually check the data, leave empty or
%                 set to 1. To not check data, set to 0.
%

% 
% OUTPUTS:
%     SEIR     -   Information structure, updated.
%
%
%   This function will save the fitted images in a new directory called
%   "fitT1_GS" within the SEIRdir. 
%
% WEB RESOURCES:
%   web('http://white.stanford.edu/newlm/index.php/Quantitative_Imaging','-browser');
% 
% 
% EXAMPLE USAGE:
%   SEIRdir = '/biac2/wandell2/data/WMDevo/adult/109_AL/QuantitativeImaging/20110622_0582/SEIR_epi_1'
%   mrQ_fitSEIR_T1_ver2(SEIRdir,[],[]);
%   
% 
% Written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University
%  The function was edited by Aviv Mezer since 2011 this verstion is from
%  2016
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016


%% Check Inputs

% This makes many assumptions about the paths. This should be flexible. All
% we need is the data - the input should just point to that. The function
% can then move one directory up from there and create the fit* directory.

if notDefined('SEIRdir') || ~exist(SEIRdir,'dir')
    SEIRdir = uigetdir(pwd,'Select your SEIR base data directory');
end

seirDataDir = fullfile(SEIRdir,'data');

if notDefined('outDir')
    outDir = fullfile(SEIRdir,'fitT1_GS');
end

if ~exist(outDir,'dir'), mkdir(outDir); end


if notDefined('checkData') || isempty(checkData) || checkData > 1
    checkData = 1;
end

% Which algorithm to use: 
%real nonlinear least squares polarity restoration. 

method = 'NLSPR'; % magnitude data

close all


%% Load the data file (this will load 'data', 'xform', and 'extra')

% We should now be pointed to this file, and not to the directory. It is
% much easier to work with filenames than with directories.

filename = 'SEIR_Dat';
loadStr = fullfile(seirDataDir, filename);

outName = ['T1Fit' method '_' filename];  
saveStr = fullfile(outDir, outName);


%% Perform T1 fits 

% Estimates T1 together with:
%   NLSPR: a and b parameters to fit the data to |a + b*exp(-TI/T1)|
T1FitExperimentData(loadStr, saveStr, method, checkData);


% Load the fitted data 
% ** We have to do it this way because T1FitExperimentData does not return
% the data as desired.
load(saveStr)
load(loadStr)


%% Save Nifti

% The 'll_T1' variable comes from the function T1FitExperimentData - it
% could be returned from there instead of loading the file containing it. 

T1file=[saveStr '_T1.nii.gz'];
resnormfile=[saveStr '_T1FitResid.nii.gz'];
dtiWriteNiftiWrapper(single(ll_T1(:,:,:,1)), xform, T1file); %#ok<NODEF>
dtiWriteNiftiWrapper(single(ll_T1(:,:,:,4)), xform, resnormfile);

SEIR.SEIR_epi_T1file=T1file;
SEIR.SEIR_epi_resnormfile=resnormfile;
SEIR.SEIR_epi_fitFile=saveStr;

%% Make a brain mask for SEIR on the M0 image
T1=ll_T1(:,:,:,1);
mask=T1>0.1 & T1<5000 & ~isnan(T1) & ~isinf(T1) ;  % where we have signal
M0= ll_T1(:,:,:,3); % this is a M0 waighted image of the fit.

M0(M0<0)=0;M0(isnan(M0))=0;M0(isinf(M0))=0; M0(M0>3*median(M0(mask)))=3*median(M0(mask));% let's clip outlier values
M0file=[saveStr '_M0.nii.gz']; %save M0 SEIR
dtiWriteNiftiWrapper(single(M0), xform, M0file); 

mmPerVox=[abs(xform(1,1)) abs(xform(2,2)) abs(xform(3,3))];
[brainMask,checkSlices] = mrAnatExtractBrain(M0, mmPerVox, 0.5); %scale strip FSL BET


out = fullfile(tempdir,'bet_tmp');
delete([out,'.img']) % remove tmp BET file
delete([out,'.hdr']) % remove tmp BET file

brainMask= brainMask & mask;
brainMask(:,:,1)=0;brainMask(:,:,end)=0;
BrainMaskfile=[saveStr '_BrainMask.nii.gz']; %save Brain mask SEIR
dtiWriteNiftiWrapper(single(brainMask), xform, BrainMaskfile); 
SEIR.SEIR_epi_M0file=M0file;
SEIR.SEIR_epi_Maskfile=BrainMaskfile;

return
