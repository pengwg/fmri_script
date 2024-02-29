%MAIN SCRIPT FOR PREPROCESSING AND CONNECTIVITY ANALYSIS
%% SET UP PARAMETERS
path_subjects = 'E:\WVU-RNI\FUS-AD-NEUROMOD\Data\';
path_spm = 'C:\spm12\spm12\'
path_conn = 'C:\conn22a\conn'
TR = 1;
FWMH = 8;
TotalReadOutTime = 36.297;
blipdir = -1;
Flag_VDM= 0; % 0 = run the step, 1 = don't run the step
Flag_Resclice_Unwarp= 0;
Flag_SliceTiming_Correction= 0;
Flag_register_functional= 0;
Flag_Segment = 0;
Flag_normalize_functional = 0;
Flag_register_greymatter = 0;
Flag_register_whitematter = 0;
Flag_register_csf = 0;
Flag_register_structural = 0;
Flag_smooth = 0;
Flag_MotionCheck= 0;
Flag_register_ROI = 0;
% end of parameter set up - Let's go! ****

%subjects list
subjects =dir(path_subjects);
subjects={subjects.name}';
subjects = setdiff(subjects,{'.';'..'});


%Session list
for i = 1:length(subjects)
    folder_of_sub = [path_subjects char(subjects(i,:)) '\'];
    cd (folder_of_sub);
    SESSIONS = dir(folder_of_sub);
    SESSIONS={SESSIONS.name}';
    SESSIONS = setdiff(SESSIONS,{'.';'..'});

%      SESSIONS = SESSIONS(2,:);

%Preprocessing function
Preprocessing_Function_FUSAD(path_subjects,subjects(i,:),SESSIONS,path_spm,TR,FWMH,blipdir, TotalReadOutTime,Flag_VDM,Flag_Resclice_Unwarp ,Flag_SliceTiming_Correction ,Flag_register_functional ,Flag_Segment,Flag_normalize_functional,Flag_register_greymatter ,Flag_register_whitematter,Flag_register_csf ,Flag_register_structural ,Flag_smooth ,Flag_MotionCheck,Flag_register_ROI )
end

%% Connectivity Script
%% SET UP PARAMETER CONNECTIVITY SCRIPT
BATCHFILENAME = ['E:\WVU-RNI\FUS-AD-NEUROMOD\Analysis\ADFUS-3Subjects_version0'];
ROIpath = 'E:\WVU-RNI\FUS-AD-NEUROMOD\ROI\Seed\';
path_greymatter = 'E:\WVU-RNI\FUS-AD-NEUROMOD\ROI\GreyMatter\rTPM_Graymatter.nii'
FunctionConnectivity_FUSAD(path_conn,path_subjects,subjects,TR,ROIpath,path_greymatter,BATCHFILENAME)

% for 90 days, subjetc 216, 218 and 221 relapsed. Missing subject 223
