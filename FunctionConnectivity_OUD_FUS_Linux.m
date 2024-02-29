function FunctionConnectivity_OUD_FUS(path_conn,path_subjects,subjects,TR,ROIpath, path_greymatter,BATCHFILENAME)
addpath path_conn



%% SETUP
%Define batch properties
clear BATCH;
BATCH.filename = BATCHFILENAME;
BATCH.Setup.isnew=1;
BATCH.Setup.RT = TR;
BATCH.Setup.nsubjects = size(subjects,1);
BATCH.Setup.acquisitiontype = 1; % continous (default value)
BATCH.Setup.analyses = [1,2,3]; % ROI to ROI and seed to voxel
BATCH.Setup.voxelmask = 2; %1-fixed template mask; 2-implicit uses subject specific mask
BATCH.Setup.voxelresolution = 3; % functional space used
BATCH.Setup.outputfiles = [0,0,1,1,1]; % write nii volumes for r-maps, p-maps and FDR-p-maps 
BATCH.Setup.conditions.missingdata   = 1; % allow subjects with missing condition data 
BATCH.Setup.rois.mask = 1;
ROIs  = dir( ROIpath ) ;
ROIs={ROIs.name}';
ROIs = setdiff(ROIs,{'.';'..'});
 for h=1:size(ROIs,1)
    BATCH.Setup.rois.files{h}= fullfile(ROIpath,ROIs(h,:)); 
 
 end
for i =1:length(subjects)
% specify targeted ROIs

SESSION = dir([path_subjects,subjects{i},filesep,'ses*']);
SESSION={SESSION.name}';
SESSION = setdiff(SESSION,{'.';'..'});
for j = 1:size(SESSION,1)


%   ;



%SET FOR ENCODING      

FUNC_folder = [path_subjects,subjects{i},filesep,char(SESSION(j,:)),filesep,'func', filesep];
ANAT_folder =[path_subjects,subjects{i},filesep,char(SESSION(j,:)),filesep,'anat', filesep];

RS_FILES = dir([path_subjects,subjects{i},filesep,char(SESSION(j,:)),filesep,'func' '/s8wuc*.nii']);
RS_FILES={RS_FILES.name}';
RS_FILES = setdiff(RS_FILES,{'.';'..'});
RS_FILE = [FUNC_folder, char(RS_FILES)];
   t1_list = spm_select('FPList',ANAT_folder,'^wsub*');
    t1_file = char(t1_list(contains(string(t1_list),'.nii'),:));
%     c1_file = spm_select('FPList',ANAT_folder,'^wc1');
    c2_file = spm_select('FPList',ANAT_folder,'^wc2');
    c3_file = spm_select('FPList',ANAT_folder,'^wc3');

% specificy structural data per subject
    BATCH.Setup.structurals{i}{j} = t1_file;
    BATCH.Setup.masks.Grey.files = path_greymatter ;
    BATCH.Setup.masks.Grey.dimensions=1;
    BATCH.Setup.masks.White.files{i}{j} = c2_file;
    BATCH.Setup.masks.White.dimensions=1;
    BATCH.Setup.masks.CSF.files{i}{j} = c3_file;
    BATCH.Setup.masks.CSF.dimensions=1; 

% specify conditions per subject {ncond}{nsub}{nsess}

BATCH.Setup.conditions.names = {'RS_Baseline','RS_7Days','RS_30Days','RS_90Days'}; %To adaptaccordingly in function of the number of session
 

%{nsub}{nses} 

 BATCH.Setup.functionals{i}{j} = spm_file(cellstr(spm_select('expand', RS_FILE)));%[ repmat(Nifti_outputdir_func_1,size(ls([Nifti_outputdir_func_1 ,'*moco*.nii']),1),1) ls([Nifti_outputdir_func_1 ,'*moco*.nii'])];
 BATCH.Setup.conditions.onsets{j}{i}{j} = 9; % we remove the first 9 scans to avoid any arctefact effect of BOLD signal peak due to the entrance in the MRI
 BATCH.Setup.conditions.durations{j}{i}{j} = Inf;
%  if j == 2
 BATCH.Setup.conditions.onsets{1}{i}{2} = []; 
 BATCH.Setup.conditions.durations{1}{i}{2} = [];
  BATCH.Setup.conditions.onsets{2}{i}{1} = []; 
 BATCH.Setup.conditions.durations{2}{i}{1} = [];
%  elseif j == 3

 BATCH.Setup.conditions.onsets{3}{i}{2} = []; 
 BATCH.Setup.conditions.durations{3}{i}{2} = [];
  BATCH.Setup.conditions.onsets{3}{i}{1} = []; 
 BATCH.Setup.conditions.durations{3}{i}{1} = [];
 BATCH.Setup.conditions.onsets{1}{i}{3} = []; 
 BATCH.Setup.conditions.durations{1}{i}{3} = [];
  BATCH.Setup.conditions.onsets{2}{i}{3} = []; 
 BATCH.Setup.conditions.durations{2}{i}{3} = [];
 % if j = 4
 
 BATCH.Setup.conditions.onsets{4}{i}{1} = []; 
 BATCH.Setup.conditions.durations{4}{i}{1} = [];
  BATCH.Setup.conditions.onsets{4}{i}{2} = []; 
 BATCH.Setup.conditions.durations{4}{i}{2} = [];
 BATCH.Setup.conditions.onsets{1}{i}{4} = []; 
 BATCH.Setup.conditions.durations{1}{i}{4} = [];
  BATCH.Setup.conditions.onsets{2}{i}{4} = []; 
 BATCH.Setup.conditions.durations{2}{i}{4} = [];
   BATCH.Setup.conditions.onsets{4}{i}{3} = []; 
 BATCH.Setup.conditions.durations{4}{i}{3} = [];
 BATCH.Setup.conditions.onsets{3}{i}{4} = []; 
 BATCH.Setup.conditions.durations{3}{i}{4} = [];
%  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DENOISING
%s{ncovariate}{nsub}{nses} 
BATCH.Setup.covariates.names = {'motion'};
%     % specify covariates (nuisance regressors) per subject {ncond}{nsub}{nsess}
    Nifti_outputdir = [path_subjects,subjects{i},filesep ,sprintf('%s/',string(SESSION(j,:))),filesep];
    BATCH.Setup.covariates.files{1}{i}{j} = spm_select('FPList',FUNC_folder,'rp_c.*/.txt$');% mvt
    end
end
    BATCH.Setup.done = 1; % run the SETUP step
    BATCH.Setup.overwrite = 0;
 
BATCH.Preprocessing.filter               = [0.008, 0.09]; % set the bandpass filter limits
BATCH.Preprocessing.confounds.names      = {'White Matter', 'CSF', 'motion'};
BATCH.Preprocessing.confounds.dimensions = {'1','1','6',};
BATCH.Preprocessing.confounds.deriv      = {'0','0','1'};
BATCH.Preprocessing.done                 = 1; 
 conn_batch(BATCH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST-LEVEL ANALYSIS
% the rest is done from the gui. See the resulting
% FUS_OUD_FSSEGMENTATION_slicetime.mat file
% example of a typical first level analysis, but here since we use the
% segmentation file from free surfer it is easier to do it with the gui




% seed to whole brain accumbens
BATCH.Setup.done = 0; 
BATCH.Preprocessing.done                 = 0; 
BATCH.Analysis.analysis_number = 1;
BATCH.Analysis.name = 'SeedNAC-WholebRain';
BATCH.Analysis.type = 2; % do seed to voxel analysis
BATCH.Analysis.measure = 1; % measure bivariate correlation
BATCH.Analysis.weight = 2; 
BATCH.Analysis.sources.names = ROIs;
BATCH.Analysis.sources.dimensions = repmat({1},1,length(cellstr(ROIs)));
BATCH.Analysis.sources.deriv = repmat({0},1,length(cellstr(ROIs)));
BATCH.Analysis.done = 1; % run the ANALYSIS step for ROI to ROI and seed to voxel
BATCH.vvAnalysis.done = 0; % run the ANALYSIS step for voxel to voxel
conn_batch(BATCH);


end