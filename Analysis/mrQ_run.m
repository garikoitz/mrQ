function mrQ_run(dir,outDir,inputData_spgr,inputData_seir,B1file, varArgIn)
% mrQ_run_Ver2(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir,B1file)
%  mrQ_run_Ver2(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir,B1file)
%
% This is an improved version of: mrQ_runNIMS(dir,Callproclus,refFile,outDir)
%
%    INPUT:
%
%                dir:   Directory where the nifti from NIMS are located.
%             outDir:   Directory to which the output will be saved. 
%                           (default: pwd/mrQ)
%     inputData_spgr:   The SPGR data
%     inputData_seir:   The SEIR data
%             B1file:   If empty (default), the function will calculate a 
%                          B1 file from the data SEIR and SPGR data.
%                          Alternatively, the B1 inhomogeneity map can be
%                          provided as a NIfTI. The file has to be
%                          registered to the SPGR reference file. Note that
%                          if no reference file is provided, AC-PC
%                          alignment will be defined below inside the
%                          function mrQ_initSPGR_ver2.m.
%           varargin:   every parameter that you would like to be 
%                          different from default can be given here. This
%                          input will be expected to come as a cell array
%                          and in pairs, according to the options given in
%                          the mrQ_set function. For example, if a
%                          reference file is given, the input should be
%                          {'ref', refFile} where refFile is the path to a
%                          reference image (nii.gz). Another example is the
%                          choice of using sungrid, in which case the input
%                          should be : {'sungrid',useSunGrid} where
%                          useSunGrid is 1/0 depending on whether or the
%                          user would like to use sungris (1) or not (0).
%                          the default is 0. One could also ask for both
%                          option: {'ref', refFile,'sungrid',useSunGrid}
%example= mrQ_run_Ver2(dir,outdir,[],[],[],{'lsq',1})
%   OUTPUT:
%
%       This function creates and saves the mrQ strucure to the subject's
%       directory. New directories will be created, including directories
%       for data and quantitative fits. Images will be registered to each
%       other. 
%            *  SEIR-EPI T1 will be computed (low resolution). 
%            *  SPGR T1, M0, B1 maps, and a synthetic T1-weighted image 
%                   will be computed.
%            *  T1-weighted and quantitative T1 images will be combined to 
%                   segment the brain tissue.
%`           *  PD and coil gain will be fitted from the M0 image.
%            *  Biophysical models will be applied to calculate VIP and 
%                   SIR maps.
% 
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%  2015
%
%  VAmos a meter a mano los TE y TR hasta que no arreglen esto...
% dataDir is the location where the NIfTI are saved.





% inDir = '/sni-storage/wandell/users/glerma/Documents/mrQPRUEBA/PILOT05/nifti';
% outDir = '/sni-storage/wandell/users/glerma/Documents/mrQPRUEBA/PILOT05/mrQAnalysis'; 
% inDir = '/Users/gari/gDrive/BCBL/PROYECTOS/MINI/mrQ2_MBP/nifti';
% outDir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/mrQ/PILOT05'; 
% inDir = '/bcbl/home/public/Gari/MINI/DATA/dicoms/S001_DAY2_MINITWICE_4189/qMRI/nifti';
% outDir = '/bcbl/home/public/Gari/MINI/ANALYSIS/mrQ/S001_DAY2_MINITWICE_4189'; 
% 
% inputData_spgr = struct;
% inputData_spgr.rawDir = inDir;
% %
% % A list of NIfTI names. (A unique string from the names is enough)
% inputData_spgr.name = {'FA4' 'FA10' 'FA20' 'FA30'};
% %
% % The TR of each NIfTI in the list (in msec)
% inputData_spgr.TR = [12 12 12 12];
% %
% % The TE of each NIfTI in the list (in msec)
% inputData_spgr.TE = [2.27 2.27 2.27 2.27];
% %
% % The flip angle of each NIfTI in the list (in degrees)
% inputData_spgr.flipAngle = [4 10 20 30];
% %
% % The field strength of each NIfTI in the list (in Teslas)
% inputData_spgr.fieldStrength = [3 3 3 3];
% % Next, create a structure called "inputData_seir". Set the required SEIR parameters:
% 
% %        B. Define the SEIR header info:
% %
% inputData_seir = struct;
% % dataDir is the location where the NIfTI are saved
% inputData_seir.rawDir = inDir;
% %
% % A list of NIfTI names.  (A unique string from the names is enough)
% inputData_seir.name = {'IT50'  'IT400'  'IT1200'  'IT2400'};
% %
% % The TR of each NIfTI in the list (in msec)
% inputData_seir.TR = [3000 3000 3000 3000];
% %
% % The TE of each NIfTI in the list (in msec)
% inputData_seir.TE = [49 49 49 49];
% %
% % The inversion time of each NIfTI in the list (in msec)
% inputData_seir.IT = [50 400 1200 2400];
% 
% % RUN IT
% mrQ_run(inDir, outDir, inputData_spgr, inputData_seir,[], {'sungrid', 0})

% mrQ_run(inDir, outDir, [], [],[], {'sungrid', 0})





%% I. Create the initial structure
 
if notDefined('outDir') 
     outDir = fullfile(dir,'mrQ');
end

% Creates the name of the output directory
if ~exist(outDir,'dir'); mkdir(outDir); end

% Creates the mrQ structure
mrQ = mrQ_Create(dir,[],outDir);

%     Create a file containing mrQpath, named after its 'ID' (taken from
%     its tempname). This allows for an easy use of SunGrid.
mrQ_createIDfile(mrQ);

% Set other parameters, such as SUNGRID and fieldstrength

if ~notDefined('varArgIn')
    if ~isempty(varArgIn)
        for ii = 1:2:numel(varArgIn)-1
            % Check to make sure that the argument is formatted properly
            mrQ = mrQ_Set(mrQ, varArgIn{ii}, varArgIn{ii+1});
        end
    end
end

%% II. Arrange the data
% A specific arrange function for nimsfs, nifti, or using input for user

if ~isfield(mrQ,'Arrange_Date');
    
    if (~notDefined('inputData_spgr') &&  ~notDefined('inputData_seir'))
        mrQ = mrQ_arrangeData_nimsfs(mrQ,inputData_spgr,inputData_seir);
    else
        mrQ = mrQ_arrangeData_nimsfs(mrQ);
        
    end
else
    fprintf('Data was already arranged on %s \n',mrQ.Arrange_Date)
end


%% IV. Initiate and align SPGR
%  parameters for aligning SPGR

%load(name);

if isfield(mrQ,'SPGR_init_done');
else
    mrQ.SPGR_init_done=0;
end

if     mrQ.SPGR_init_done==0
    
    % Keeps track of the variables we use.
    % For details, look inside the function.
    [~, ~, ~,~,~, mrQ]=mrQ_initSPGR(mrQ.SPGR,mrQ.refIm,mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
    mrQ.SPGR_init_done=1;
    
    save(mrQ.name,'mrQ');
    fprintf('\n  init SPGR - done!           \n');
else
    fprintf(' \n Loading init SPGR data            \n');
    
end

%%  V. Fit SPGR PD

if ~isfield(mrQ,'SPGR_LinearT1fit_done');
    
    mrQ.SPGR_LinearT1fit_done=0;
end

% clobber is implemented inside (we can add this to the inputs)
if (mrQ.SPGR_LinearT1fit_done==0);
    
    [mrQ]=mrQfit_T1M0_Lin(mrQ);
    
    mrQ.SPGR_LinearT1fit_done=1;
    
    save(mrQ.name,'mrQ');
    
    fprintf('\n Fit linear T1 SPGR  - done!              \n');
else
    fprintf('\n Loading linearly fitted SPGR T1                \n');
    
end

%% III. Perform SEIR fit

if notDefined('B1file')
    % Checks if B1 was defined by the user.
    % If not, we will use the SEIR data to map it.
    
    if isfield(mrQ,'SEIR_done');
    else
        mrQ.SEIR_done=0;
    end
    
    if (mrQ.SEIR_done==0);
        mrQ=mrQ_SEIR(mrQ);
        
    else
        fprintf('\n Loading previously fitted SEIR data \n');        
    end
end
%% VII. Build B1

if notDefined('B1file')
    % Checks if B1 was defined by the user.
    
    if ~isfield(mrQ,'B1Build_done');
        mrQ.B1Build_done=0;
    end
    
    if ( mrQ.B1Build_done==0)
        
        
        mrQ=mrQ_B1_LR(mrQ);
        
        mrQ.B1Build_done=1;
        save(mrQ.name,'mrQ');
        fprintf('\n Building B1 - done!       \n');
        
    else
        fprintf(['Using existing B1 map file '   mrQ.B1FileName        '  \n']);
        
    end
else
    mrQ.B1FileName=B1file;
end

%% VIII. T1M0 fit with B1

if ~isfield(mrQ,'SPGR_T1fit_done');
    mrQ.SPGR_T1fit_done=0;
end

if ( mrQ.SPGR_T1fit_done==0)
    
    
    mrQ=mrQ_T1M0_Fit(mrQ);
    mrQ.SPGR_T1fit_done=true;
    save(mrQ.name,'mrQ');
    
    fprintf('\n fit T1 SPGR  - done!              \n');
    
else
    
    fprintf('\n Using previously fitted SPGR T1                \n');
    
end

%%  Segmentation needed for PD fit
% Prefer to PD fit 
% 1. Get a segmentation (need freesurfer output) 
% 2. Get CSF 
% 3. Make a M0 file for the coils

%% IX. Create the synthetic T1-weighted images and save them to disk

if ~isfield(mrQ,'synthesis')
    mrQ.synthesis=0;
end
if mrQ.synthesis==0
    
    [mrQ.SegInfo.T1wSynthesis_MOT1,mrQ.SegInfo.T1wSynthesis_T1] =mrQ_T1wSynthesis1(mrQ,[],[],mrQ.HeadMask);
    
    mrQ.synthesis=1;
    save(mrQ.name, 'mrQ')
    
    fprintf('\n Synthesis of T1  - done!              \n');
else
    fprintf('\n Using previously synthesized T1              \n');
end


%% X. Segmentation and CSF

if ~isfield(mrQ,'segmentation');
    mrQ.segmentation=0;
end

if mrQ.segmentation==0;
    
    %     default- FSL segmentation
    if (mrQ.runfreesurfer==0 && ~isfield(mrQ,'freesurfer'))
        % Segment the T1w by FSL (step 1) and get the tissue mask (CSF-WM-GM) (step 2)
        %         mrQ=mrQ_segmentT1w2tissue(mrQ);
        mrQ=mrQ_Seg_kmeans(mrQ,mrQ.BrainMask);
        mrQ.segmentation=1;
        
        %      run FreeSurfer: it is slow and needs extra definitions.
    elseif (mrQ.runfreesurfer==1)
        mrQ=mrQ_Complitfreesurfer(mrQ);
        mrQ.segmentation=1;
        
        %      use an uploaded freesurfer nii.gz
    elseif   isfield(mrQ,'freesurfer');
        [mrQ.SegInfo]=mrQ_CSF(mrQ.spgr_initDir,mrQ.freesurfer,[],mrQ.AnalysisInfo);
        mrQ.segmentation=1;       
    end
    
    save(mrQ.name,'mrQ');
    fprintf('\n Segmentation and CSF  - done!              \n');
else
    fprintf('\n Using previously segmented data              \n');
    
    
end

%% XI. Fitting PD from M0

if ~isfield(mrQ,'PDdone')
    mrQ.PDdone=0;
end
if mrQ.PDdone==0
    
   mrQ=mrQ_M0_ToPD(mrQ);
    mrQ.PDdone=1;
    save(mrQ.name, 'mrQ')
    fprintf('\n Calculation of PD from M0  - done!              \n');
else
    fprintf('\n Using previously calculated PD              \n');
end


%% XII. Calculate VIP, TV,  SIR and synthetic T1w 

if ~isfield(mrQ,'VIP_WF_done')
    mrQ.VIP_WF_done=0;
end

% if (mrQ.VIP_WF_done==0)
%     fprintf('\n Calculate VIP, TV and SIR form T1 and WF maps               \n');
%     
%     [mrQ] = mrQ_WF(mrQ);
%     
%     % GLU: ERROR: Field assignment to a non-structure array object.
%     % [mrQ.AnalysisInfo, mrQ] = mrQ_VIP(mrQ);
%     % I edited mrQ_VIP.m too
%     % GLU end.
%     mrQ = mrQ_VIP(mrQ);
%     
%     mrQ.VIP_WF_done=1;
%     save(mrQ.name,'mrQ');
%     
%     fprintf('\n Calculation of VIP, MTV and SIR  - done!              \n');
%     %
%     % XIII. Create a series of synthetic T1w images
% 
%     [mrQ.T1w_file,mrQ.T1w_file1] =mrQ_T1wSynthesis1(mrQ);
%     fprintf('\n Calculation synthetic T1w- done!              \n');
% 
% else 
%      fprintf('\n Using previously calculated VIP, MTV, SIR and synthetic T1w             \n');
% end

if (mrQ.VIP_WF_done==0)
    fprintf('\n Calculate VIP, TV and SIR form T1 and WF maps               \n');
    
    [mrQ] = mrQ_WF(mrQ);
    
    
    [mrQ.AnalysisInfo, mrQ] = mrQ_VIP(mrQ);
    
    mrQ.VIP_WF_done=1;
    save(mrQ.name,'mrQ');
    
    fprintf('\n Calculation of VIP, MTV and SIR  - done!              \n');
    %
    % XIII. Create a series of synthetic T1w images

    [mrQ.T1w_file,mrQ.T1w_file1] =mrQ_T1wSynthesis1(mrQ);
    fprintf('\n Calculation synthetic T1w- done!              \n');

else 
     fprintf('\n Using previously calculated VIP, MTV, SIR and synthetic T1w             \n');
end


%% XIV. Organize the OutPut directory
mrQ=mrQ_arrangeOutPutDir(mrQ);
mrQ_deleteIDfile(mrQ);% delete the temporary ID file stored in mrQ/sge_subjects

%done
mrQ.AnalysisDone=1;
mrQ.AnalysisDoneDate=date;
%save
save(mrQ.name,'mrQ');
fprintf('\n done !!! \n')
