clear all
close all
% %Add SPM path
addpath 'C:\spm12\spm12' %To adapt accordingly
%Add CONN toolbox path
addpath 'C:\conn20b\conn'%To adapt accordingly

%Load subject folders
path_subjects='C:\Users\natha\Downloads\FUS_transferToNIDA(1)\FUS_transferToNIDA\';%To adapt accordingly
Subjects=dir(path_subjects);
Subjects={Subjects.name}';
Subjects = setdiff(Subjects,{'.';'..'});
Subjects = Subjects(end-9:end,:)

RT = 1;


%% SETUP
%Define batch properties
clear BATCH;
BATCH.filename = ['D:\WVU-RNI\FUS-OUD-resting_state\Analysis\FUS_FINALNOVEMBER_complete'];
BATCH.Setup.isnew=1;
BATCH.Setup.RT = RT;
BATCH.Setup.nsubjects = size(Subjects,1);
BATCH.Setup.acquisitiontype = 1; % continous (default value)
BATCH.Setup.analyses = [1,2,3]; % ROI to ROI and seed to voxel
BATCH.Setup.voxelmask = 2; %1-fixed template mask; 2-implicit uses subject specific mask
BATCH.Setup.voxelresolution = 3; % functional space used
BATCH.Setup.outputfiles = [0,0,1,1,1]; % write nii volumes for r-maps, p-maps and FDR-p-maps 
BATCH.Setup.conditions.missingdata   = 1; % allow subjects with missing condition data 
BATCH.Setup.rois.mask = 1;
for i =1:length(Subjects)
% specify targeted ROIs
SESSION = dir([path_subjects,Subjects{i},filesep,'ses*']);
SESSION={SESSION.name}';

for j = 1:size(SESSION,1)

ROIpath = [path_subjects char(Subjects(i,:))  '\' char(SESSION(j,:)) '\ROI\ROI_FCANALYSIS_MNI\']; % si far don't have the segmentation for everyone so take the segmentation of the second session for everything
ROIs  = ls([path_subjects char(Subjects(i,:))  '\' char(SESSION(j,:)) '\ROI\ROI_FCANALYSIS_MNI\w*.nii']) ;
%   ;


 for h=1:size(ROIs,1)
    BATCH.Setup.rois.files{h}{i}{j}= fullfile(ROIpath,ROIs(h,:)); 
 
%  if i == length(Subjects) && j == size(SESSION,1)
% 
%      INDEX = strfind(ROIs(h,:),'_');
%      INDEX = INDEX(2);
%      BATCH.Setup.rois.names{h} = cellstr(ROIs(h,1:INDEX-1));
%  end
 end

%SET FOR ENCODING      

FUNC_folder = [path_subjects,Subjects{i},filesep,char(SESSION(j,:)),filesep,'func', filesep];
ANAT_folder =[path_subjects,Subjects{i},filesep,char(SESSION(j,:)),filesep,'anat', filesep];
RS_FILE = [FUNC_folder, ls([path_subjects,Subjects{i},filesep,char(SESSION(j,:)),filesep,'func' '\swuc*.nii'])];
   t1_list = spm_select('FPList',ANAT_folder,'^wsub*');
    t1_file = char(t1_list(contains(string(t1_list),'.nii'),:));
    c1_file = spm_select('FPList',ANAT_folder,'^wc1');
    c2_file = spm_select('FPList',ANAT_folder,'^wc2');
    c3_file = spm_select('FPList',ANAT_folder,'^wc3');

% specificy structural data per subject
    BATCH.Setup.structurals{i}{j} = t1_file;
    BATCH.Setup.masks.Grey.files{i}{j} = c1_file;
    BATCH.Setup.masks.Grey.dimensions=1;
    BATCH.Setup.masks.White.files{i}{j} = c2_file;
    BATCH.Setup.masks.White.dimensions=1;
    BATCH.Setup.masks.CSF.files{i}{j} = c3_file;
    BATCH.Setup.masks.CSF.dimensions=1; 

% specify conditions per subject {ncond}{nsub}{nsess}

BATCH.Setup.conditions.names = {'RS_Baseline','RS_7Days','RS_30Days'};
 

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
%  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DENOISING
%s{ncovariate}{nsub}{nses} 
BATCH.Setup.covariates.names = {'motion'};
%     % specify covariates (nuisance regressors) per subject {ncond}{nsub}{nsess}
    Nifti_outputdir = [path_subjects,Subjects{i},filesep ,sprintf('%s\',string(SESSION(j,:))),filesep];
    BATCH.Setup.covariates.files{1}{i}{j} = spm_select('FPList',FUNC_folder,'rp_c.*\.txt$');% mvt
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


brainRegions={'caudalanteriorcingulate', 'rostralanteriorcingulate',  'posteriorcingulate' ,  'lateralorbitofrontal', ...
'medialorbitofrontal' , 'caudalmiddlefrontal'  ,  'rostralmiddlefrontal' ,'frontalpole' , 'insula' , 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Amygdala', 'Accumbens', 'parsorbitalis', 'parstriangularis', 'parsopercularis'};

brainRegions_Cortical = {'caudalanteriorcingulate', 'rostralanteriorcingulate',  'posteriorcingulate' ,  'lateralorbitofrontal', ...
'medialorbitofrontal' , 'caudalmiddlefrontal'  ,  'rostralmiddlefrontal' ,'frontalpole' , 'insula' , 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Amygdala', 'Accumbens','parsorbitalis', 'parstriangularis', 'parsopercularis','pericalcarine'};

brainRegions_Subcortical = { 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Amygdala', 'Accumbens'};
Accumbens = {'Accumbens'};

% seed to whole brain accumbens
BATCH.Setup.done = 0; 
BATCH.Preprocessing.done                 = 0; 
BATCH.Analysis.analysis_number = 1;
BATCH.Analysis.name = 'SeedNAC-WholebRain';
BATCH.Analysis.type = 2; % do seed to voxel analysis
BATCH.Analysis.measure = 1; % measure bivariate correlation
BATCH.Analysis.weight = 2; 
BATCH.Analysis.sources.names = ROIs(ismember(ROIs,Accumbens ));
BATCH.Analysis.sources.dimensions = repmat({1},1,length(cellstr(ROIs(ismember(ROIs,Accumbens )))));
BATCH.Analysis.sources.deriv = repmat({0},1,length(cellstr(ROIs(ismember(ROIs,Accumbens )))));
BATCH.Analysis.done = 1; % run the ANALYSIS step for ROI to ROI and seed to voxel
BATCH.vvAnalysis.done = 0; % run the ANALYSIS step for voxel to voxel
conn_batch(BATCH);

%Cortical
BATCH.Setup.done = 0; 
BATCH.Preprocessing.done                 = 0; 
BATCH.Analysis.analysis_number = 1;
BATCH.Analysis.name = 'ROI_cortical';
BATCH.Analysis.type = 2; % do seed to voxel analysis
BATCH.Analysis.measure = 1; % measure bivariate correlation
BATCH.Analysis.weight = 2; 
 BATCH.Analysis.sources.names = ROIs(ismember(ROIs,brainRegions_Cortical ));
BATCH.Analysis.sources.dimensions = repmat({1},1,length(cellstr(ROIs(ismember(ROIs,brainRegions_Cortical )))));
BATCH.Analysis.sources.deriv = repmat({0},1,length(cellstr(ROIs(ismember(ROIs,brainRegions_Cortical )))));
BATCH.Analysis.done = 1; % run the ANALYSIS step for ROI to ROI and seed to voxel
BATCH.vvAnalysis.done = 0; % run the ANALYSIS step for voxel to voxel
conn_batch(BATCH);


%SubCortical
BATCH.Setup.done = 0; 
BATCH.Preprocessing.done                 = 0; 
BATCH.Analysis.analysis_number = 1;
BATCH.Analysis.name = 'ROI_Subcortical';
BATCH.Analysis.type = 2; % do seed to voxel analysis
BATCH.Analysis.measure = 1; % measure bivariate correlation
BATCH.Analysis.weight = 2; 
 BATCH.Analysis.sources.names = ROIs(ismember(ROIs,brainRegions_Subcortical ));
BATCH.Analysis.sources.dimensions = repmat({1},1,length(cellstr(ROIs(ismember(ROIs,brainRegions_Subcortical )))));
BATCH.Analysis.sources.deriv = repmat({0},1,length(cellstr(ROIs(ismember(ROIs,brainRegions_Subcortical )))));
BATCH.Analysis.done = 1; % run the ANALYSIS step for ROI to ROI and seed to voxel
BATCH.vvAnalysis.done = 0; % run the ANALYSIS step for voxel to voxel
conn_batch(BATCH);

