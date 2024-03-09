% MAIN SCRIPT FOR PREPROCESSING AND CONNECTIVITY ANALYSIS
clear

% TO ADAPT ACCORDINGLY
path_subjects = '~/Work/fusOUD/fmri/';
path_subjects = '/tm/WVU-RNI/FUS-OUD/PreprocessedData/';

% MR parameters
TR = 1;
FWMH = 4;
TotalReadOutTime = 36.297;
blipdir = -1;

SPM_directives.do_VDM = true;
SPM_directives.do_Resclice_Unwarp = true;
SPM_directives.do_SliceTiming_Correction = true;
SPM_directives.do_register_functional = true;
SPM_directives.do_Segment = true;
SPM_directives.do_normalize_functional = true;
SPM_directives.do_register_greymatter = true; 
SPM_directives.do_register_whitematter = true;
SPM_directives.do_register_csf = true;
SPM_directives.do_register_structural = true;
SPM_directives.do_smooth = true;
SPM_directives.do_MotionCheck = true;
SPM_directives.do_register_ROI = true;

subjects = dir([path_subjects '/sub*']);
% TO ADAPT ACCORDINGLY (for selective processing)
% subjects = subjects(end-8:end, :); 

% ROI conversion from MGH to NII (to extract NAC only)
% conn_mgh2nii(path and name of file to convert here)

for i = length(subjects)
    SESSIONS = dir([path_subjects subjects(i).name '/ses-*']);

    spm_preprocess(subjects(i).name, SESSIONS, TR, FWMH, blipdir, TotalReadOutTime, SPM_directives)
end

% Connectivity Script
% subjects = subjects(~contains (subjects, '220'),:);

BATCHFILENAME = [path_subjects 'Analysis/FUS_FINALtrial_90Days2'];
ROIpath = 'ROI/Seed/';
path_greymatter = 'ROI/GreyMatter/rTPM_Graymatter.nii';
conn_FUSOUD(path_subjects, {subjects.name}, TR, ROIpath, path_greymatter, BATCHFILENAME)

% for 90 days, subjetc 216, 218 and 221 relapsed. Missing subject 223
