function mrQ=mrQ_SEIR(mrQ)
% we found the registration is often failing (causing a bias in the T1 fit
% in following steps). A possible solution is to re-align (and then re-fit)
% the SEIR, with a different IT image as the target of the registration. To
% monitor that, we register according to the first TI, and then check the
% registration using the funcction mrQ_QuantAnts. The check is not optimal
% and is using a somewhat arbitrary threshold. if the check doesn't pass
% the threshold - we reallign, fit and register spgr to the new fit - and
% then check again. if non of the trials pass our threshold, we ask the
% user to mannually check the registration, and if it is satisfying one can
% easilly bypass out checks using, for example, rerun mrQ with the
% additional input {'seir_done',1}: 
% mrQ_run(inputDir,outDir,[],[],[], {'seir_done',1}),
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel, Copyright 2016

if isfield(mrQ,'AntSQantTresh');  
else
    mrQ.AntSQantTresh=0.65;
end

AntsTmpPath = fullfile(mrQ.spgr_initDir, 'tmpAnts');
SEIRFitTmpPath = fullfile(mrQ.SEIRepiDir, 'tmpSEIRFit');
SEIR_fit_dir = fullfile(mrQ.SEIRepiDir, 'fitT1_GS');
if ~exist(AntsTmpPath,'dir'), mkdir(AntsTmpPath); end 
if ~exist(SEIRFitTmpPath,'dir'), mkdir(SEIRFitTmpPath); end 
if ~exist(SEIR_fit_dir, 'dir'), mkdir(SEIR_fit_dir); end 

mrQ.AntSQantToHigh=inf;   

% <<<<<<< HEAD: GLU: just in case maintaining the conflicts as comments
for ii=1: length(mrQ.inputdata_seir.IT) %% in each iteration, align the SEIR images with a different image as the target image
    % Keeps track of the variables we use.
    % For details, see inside the function.
    [~, ~, ~, mrQ.SEIRsaveData]=mrQ_initSEIR(mrQ,mrQ.SEIRepiDir,mrQ.alignFlag,ii);
    
    SEIR_epi_AlignIm=mrQ.inputdata_seir.name{ii};  %keep the image that is the target image in the current iteration. 
    
    [mrQ]=mrQ_fitSEIR_T1(mrQ.SEIRepiDir,SEIRFitTmpPath,0,mrQ);
    
    save(mrQ.name,'mrQ');
    fprintf('Fit SEIR  - done! \n');
    

   
    %%   Register high-resolution EPI image to low-resolution aligned T1 image
    
    %mrQ_NLANTS_warp_SPGR2EPI_RB(AnalysisInfo,SET1file,t1fileHM,flipAngles,outDir,AlignFile)
    %
    % if ~isfield(mrQ,'SPGR_EPI_align_done');
    %     mrQ.SPGR_EPI_align_done=0;
    % end
    %
    % if ( mrQ.SPGR_EPI_align_done==0)
    
    % create a temp folder for the ...

      
    if isfield(mrQ,'antsmet')
        if mrQ.antsmet=='1'
            % fitting based on the brain only- no head
            [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIR_epi_T1file,mrQ.T1_LFit,mrQ.SEIR_epi_Maskfile,AntsTmpPath,{mrQ.T1_LFit_HM});
        elseif mrQ.antsmet=='2'
            % fitting without a mask on the epi
            [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIR_epi_T1file,mrQ.T1_LFit_HM,[],AntsTmpPath,{mrQ.T1_LFit_HM});
        end
    else
        [WARP_SPGR_EPI,  T1_spgr_epi]= mrQ_NLANTS_SPGR2EPI(mrQ.SEIR_epi_T1file,mrQ.T1_LFit_HM,mrQ.SEIR_epi_Maskfile,AntsTmpPath,{mrQ.T1_LFit_HM});
        
    end
   
      
    MovingScaleConstat=1000;% to compare SPGR T1 in sec to SEIR T1 in ms.
    % check the registration:
    T1_spgr_epi = T1_spgr_epi{1};
    [mrQ.AntSQant]=mrQ_QuantAnts(mrQ.SEIR_epi_T1file,T1_spgr_epi,MovingScaleConstat);

    % We will keep the files only if the are the current best Ants registration
     if min(mrQ.AntSQantToHigh) > mrQ.AntSQant 
        
        cmd= [ 'cp ' T1_spgr_epi ' ' mrQ.spgr_initDir '/.'  ];
        system(cmd);
        cmd= [ 'cp ' WARP_SPGR_EPI '* ' mrQ.spgr_initDir '/.'  ];
        system(cmd);
        

        
        [~, Name]=fileparts(WARP_SPGR_EPI);
        mrQ.Ants_Info.WARP_SPGR_EPI = fullfile(mrQ.spgr_initDir,Name);
         
        [~, Name, ext]=fileparts(T1_spgr_epi);
        Name = [Name, ext];
        mrQ.Ants_Info.T1_spgr_epi=fullfile(mrQ.spgr_initDir,Name);
        
        % Saving the epi alignment parameters
        
        

        mrQ.SEIR_epi_AlignIm=SEIR_epi_AlignIm;
        
        cmd= [ 'cp ' mrQ.SEIR_epi_T1file ' ' SEIR_fit_dir '/.'  ];
        system(cmd);
        SEIR_epi_T1file = mrQ.SEIR_epi_T1file;
        
        cmd= [ 'cp ' mrQ.SEIR_epi_resnormfile ' ' SEIR_fit_dir '/.'  ];
        system(cmd);
        SEIR_epi_resnormfile = mrQ.SEIR_epi_resnormfile;
        
        cmd= [ 'cp ' mrQ.SEIR_epi_fitFile '.mat ' SEIR_fit_dir '/.'  ];
        system(cmd);
        SEIR_epi_fitFile = mrQ.SEIR_epi_fitFile;
        
        cmd= [ 'cp ' mrQ.SEIR_epi_M0file ' ' SEIR_fit_dir '/.'  ];
        system(cmd);
        SEIR_epi_M0file = mrQ.SEIR_epi_M0file;
        
        cmd= [ 'cp ' mrQ.SEIR_epi_Maskfile ' ' SEIR_fit_dir '/.'  ];
        system(cmd);
        SEIR_epi_Maskfile = mrQ.SEIR_epi_Maskfile;
                     
    end

    
    if mrQ.AntSQant< mrQ.AntSQantTresh;
        break
    else
        % if the registration quality did not pass the threshold, we will
        % keep the values of the registration quality and move on to the
        % next iteration. 
        
%         mrQ.Ants_Info.WARP_SPGR_EPI = best_WARP_SPGR_EPI;
%         mrQ.Ants_Info.T1_spgr_epi = best_T1_spgr_epi;
        mrQ.AntSQantToHigh(ii)=mrQ.AntSQant;
        mrQ.SEIR_epi_AlignImToHigh{ii}=SEIR_epi_AlignIm;
        save(mrQ.name,'mrQ');
    end
end

%saving the correct SEIR files
[~, Name, ext]=fileparts(SEIR_epi_T1file);
Name = [Name, ext];
mrQ.SEIR_epi_T1file = fullfile(SEIR_fit_dir,Name);

[~, Name, ext]=fileparts(mrQ.SEIR_epi_resnormfile);
Name = [Name, ext];
mrQ.SEIR_epi_resnormfile = fullfile(SEIR_fit_dir,Name);

[~, Name]=fileparts(SEIR_epi_fitFile);
mrQ.SEIR_epi_fitFile = fullfile(SEIR_fit_dir,Name);

[~, Name, ext]=fileparts(SEIR_epi_M0file);
Name = [Name, ext];
mrQ.SEIR_epi_M0file = fullfile(SEIR_fit_dir,Name);

[~, Name, ext]=fileparts(SEIR_epi_Maskfile);
Name = [Name, ext];
mrQ.SEIR_epi_Maskfile = fullfile(SEIR_fit_dir,Name);

% deleting temp folders
cmd = ['rm -Rf ', AntsTmpPath];
system(cmd);
cmd = ['rm -Rf ', SEIRFitTmpPath];
system(cmd);

        
if mrQ.AntSQant> mrQ.AntSQantTresh;
    
    % saving the correct AntSQant value
    mrQ.AntSQant = min(mrQ.AntSQantToHigh);
    save(mrQ.name,'mrQ');
% ======= GLU commented
% mrQ=mrQ_Call_AntsAlign_forSEIR_SPGR(mrQ);  GLU commented
% if mrQ.Ants_Info.QuantAntsScore > mrQ.QuantAntsThresh;  GLU commented
   
% >>>>>>> c8566d45979cd83ffe8e86fdcd6a2e6d89ba1f43 GLU commented
    error('we can not trust the EPI-SPGR registration \nPlease manually check the registration between \n %s and \n %s \n If it is ok, manually change mrQ.SEIR_done to be =1.', mrQ.SEIR_epi_T1file,mrQ.Ants_Info.T1_spgr_epi)
end


mrQ.SEIR_done=1;

% mrQ.SPGR_EPI_align_done=1;

save(mrQ.name,'mrQ');
fprintf('\n Alignment of EPI to T1  - done!              \n');




