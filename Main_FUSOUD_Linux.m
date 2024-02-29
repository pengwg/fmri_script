%MAIN SCRIPT FOR PREPROCESSING AND CONNECTIVITY ANALYSIS

path_subjects = 'C:/Users/natha/Downloads/FUS_transferToNIDA(1)/FUS_transferToNIDA/';
path_spm = 'C:/spm12/spm12/'
TR = 1;
FWMH = 4;
TotalReadOutTime = 36.297;
blipdir = -1;

%subjects list
subjects =dir(path_subjects);
subjects={subjects.name}';
subjects = setdiff(subjects,{'.';'..'});
subjects = subjects(end-8:end,:); %TO ADAPT ACCORDINGLY (to have the list of subjects you want to preprocess/analyze with the connectivity analysis)

% ROI conversion from MGH to NII (to extract NAC only)
%conn_mgh2nii(path and name of file to convert here)

%Session list
for i = length(subjects)
    folder_of_sub = [path_subjects char(subjects(i,:)) '/'];
    cd (folder_of_sub);
    SESSIONS = dir(folder_of_sub);
    SESSIONS={SESSIONS.name}';
    SESSIONS = setdiff(SESSIONS,{'.';'..'});

    SESSIONS = SESSIONS(4,:);

%Preprocessing function
Preprocessing_Function(path_subjects,subjects(i,:),SESSIONS,path_spm,1,4,-1, 36.297,0,0 ,0 ,0 ,0,0,0 ,0,0 ,0 ,0 ,0,0 )
end

%Connectivity Script
subjects = subjects(~contains (subjects, '220'),:);
BATCHFILENAME = ['E:/WVU-RNI/FUS-OUD-resting_state/Analysis/FUS_FINALtrial_90Days2'];
ROIpath = 'E:/WVU-RNI/FUS-OUD-resting_state/Analysis/FUS_FINALNOVEMBER_complete/ROI/ROI_UnionNAC/ROITOUSEFORSEEDBASED/';
path_greymatter = 'E:/WVU-RNI/FUS-OUD-resting_state/Analysis/FUS_FINALNOVEMBER_complete/ROI/ROI_UnionNAC/rTPM_Graymatter.nii'
FunctionConnectivity_OUD_FUS('C:/conn22a/conn',path_subjects,subjects,TR,ROIpath,path_greymatter,BATCHFILENAME)

% for 90 days, subjetc 216, 218 and 221 relapsed. Missing subject 223
