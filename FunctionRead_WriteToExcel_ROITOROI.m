% read for each subject and each session their connectivity matrix and take the connectivity
% strength between NAC and the rest to past them into an excel file 
FLAG_WRITEROI_MNI = 0;
EXCEL_FLAG = 0 ;
FLAG_POSTHOC = 1;
FLAG_CONTROL = 1;
Path_results_Cortical ='D:\WVU-RNI\FUS-OUD-resting_state\Analysis\FUS_FINALNOVEMBER_complete\results\firstlevel\ROI_cortical\'
Path_results_Subcortical ='D:\WVU-RNI\FUS-OUD-resting_state\Analysis\FUS_FINALNOVEMBER_complete\results\firstlevel\ROI_subcortical\'
Path_Results_PostHoc = 'D:\WVU-RNI\FUS-OUD-resting_state\Analysis\FUS_FINALNOVEMBER_complete\results\firstlevel\PostHoc\'
Path_results_Pericalcarine_Subcortical = 'D:\WVU-RNI\FUS-OUD-resting_state\Analysis\FUS_FINALNOVEMBER_complete\results\firstlevel\ROI_ControlSubCorticalPericalcarine\';
SESSION = {'sess-00','sess-07','sess-30'};
NAME ='Data_FS_231123.xlsx' %name of this excel file
if FLAG_WRITEROI_MNI
%% isolate NAC for each participant and putting in MNI
path_spm ='C:\spm12\spm12\';%TO ADAPT ACCORDINGLY
path_subjects ='C:\Users\natha\Downloads\FUS_transferToNIDA(1)\FUS_transferToNIDA\'; %TO ADAPT ACCORDINGLY
subjects =dir(path_subjects);
subjects={subjects.name}';
subjects = setdiff(subjects,{'.';'..'});
subjects = subjects(end-11:end,:); %TO ADAPT ACCORDINGLY (to have the list of subjects you want to preprocess)
INDEX_ROI = {'left_Accumbens','26';'right_Accumbens','58';'left_caudalanteriorcingulate','1002';'right_caudalanteriorcingulate','2002';...
 'left_rostralanteriorcingulate','1026';'right_rostralanteriorcingulate','2026';'left_posteriorcingulate','1023';'right_posteriorcingulate','2023';...
  'left_lateralorbitofrontal','1012';'right_lateralorbitofrontal','2012';'left_medialorbitofrontal','1014';'right_medialorbitofrontal','2014';...
  'left_caudalmiddlefrontal','1003';'right_caudalmiddlefrontal','2003';'left_rostralmiddlefrontal','1027';'right_rostralmiddlefrontal','2027';...
  'left_frontalpole','1032';'right_frontalpole','2032';'left_insula','1035';'right_insula','2035';'left_Thalamus','10';'right_Thalamus','49';...
  'left_Caudate','11';'right_Caudate','50';'left_Putamen','12';'right_Putamen','51';'left_Pallidum','13';'right_Pallidum','52';...
  'left_Amygdala','18';'right_Amygdala','54';'left_ifg_parsopercularis','1018';'right_ifg_parsopercularis','2018';'left_ifg_parsorbitalis','1019';'right_ifg_parsorbitalis','2019';'left_ifg_parstriangularis','1020';'right_ifg_parstriangularis','2020';...
  'left_pericalcarine','1021';'right_pericalcarine','2021';'parahippocampus','1016';'right_parahippocampus', '2016'};
for i = 1:length(subjects)
    folder_of_sub = [path_subjects char(subjects(i,:)) '\'];
    cd (folder_of_sub);
    SESSIONS = dir(folder_of_sub);
    SESSIONS={SESSIONS.name}';
    SESSIONS = setdiff(SESSIONS,{'.';'..'});

    for j = 1:size(SESSIONS,1);
path_to_ROI = [folder_of_sub char(SESSIONS(j,:) ) '\ROI\'];
path_to_anat = [folder_of_sub char(SESSIONS(j,:) ) '\anat\'];

cd(path_to_ROI);
mkdir('ROI_FCANALYSIS');mkdir('ROI_FCANALYSIS_MNI')
cd('ROI_FCANALYSIS');

for k = 1:size(INDEX_ROI,1)
V = spm_vol([path_to_ROI ls([path_to_ROI 'a*.nii'])]);
Y = spm_read_vols(V);
    s = [char(INDEX_ROI(k,1)) '_' char(subjects(i,:)) '.nii'];
    id = find(Y == str2num(char(INDEX_ROI(k,2))));
    M=zeros(V.dim);
    M(id) = 1;
    Vout = V;
    Vout.fname = s;
    Vout.dt = [2,0];
    Vout = spm_create_vol(Vout);
    spm_write_vol(Vout,M);
end

 ROI_toNormalize = ls('*.nii');
    spm('defaults','FMRI');
    spm_jobman('initcfg');
 matlabbatch{1}.spm.spatial.preproc.channel.vols = {[path_to_anat 'zipped\' ls([path_to_anat 'zipped\sub-*.nii']) ',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'C:\spm12\spm12\tpm\TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'C:\spm12\spm12\tpm\TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'C:\spm12\spm12\tpm\TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'C:\spm12\spm12\tpm\TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'C:\spm12\spm12\tpm\TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'C:\spm12\spm12\tpm\TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
matlabbatch{2}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{2}.spm.spatial.normalise.write.subj.resample = {[path_to_anat 'zipped\' ls([path_to_anat 'zipped\sub-*.nii']) ',1']};
matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run', matlabbatch);
clear matlabbatch
% 
for k = 1:size(ROI_toNormalize,1)
def_field = spm_select('FPlist', [path_to_anat 'zipped\'],'^y_sub.*\.nii$');
ROI_FILE = [pwd '\' strtrim(ROI_toNormalize(k,:))];
    spm('defaults','FMRI');
    spm_jobman('initcfg');

matlabbatch{1}.spm.spatial.normalise.write.subj.def = {def_field};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {ROI_FILE} ;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];%TO ADAPT ACCORDINGLY
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix ='w';
spm_jobman('run', matlabbatch);
clear matlabbatch
movefile([path_to_ROI 'ROI_FCANALYSIS/w' strtrim(ROI_toNormalize(k,:))],[path_to_ROI '/ROI_FCANALYSIS_MNI/w' strtrim(ROI_toNormalize(k,:))]);
    end
    end
end
end
%%
if EXCEL_FLAG
path_subjects='C:\Users\natha\Downloads\FUS_transferToNIDA(1)\FUS_transferToNIDA\';%To adapt accordingly
subjects=dir(path_subjects);
subjects={subjects.name}';
subjects = setdiff(subjects,{'.';'..'});
subjects = subjects(end-9:end,:);%To adapt accordingly to have the list of subjects to analyze

% brainRegions={'caudalanteriorcingulate','rostralanteriorcingulate', 'posteriorcingulate' , 'lateralorbitofrontal', ...
%    'medialorbitofrontal' ,'caudalmiddlefrontal'  , 'rostralmiddlefrontal' ,'frontalpole' ,'insula' ,'Thalamus','Caudate','Putamen','Pallidum','Amygdala','Accumbens'};



count = 1;
NAME_ROI = {'wleft_Accumbens_sub-223-FUS','3';'wleft_Amygdala_sub-223-FUS','2';'wleft_Caudate_sub-223-FUS','2';'wleft_Pallidum_sub-223-FUS','2';...
   'wleft_Putamen_sub-223-FUS','2';'wleft_Thalamus_sub-223-FUS','2';'wleft_caudalanteriorcingulate_sub-223-FUS','1';'wleft_caudalmiddlefrontal_sub-223-FUS','1';...
   'wleft_frontalpole_sub-223-FUS','1';'wleft_ifg_parsopercularis_sub-223-FUS','1';'wleft_ifg_parsorbitalis_sub-223-FUS','1';'wleft_ifg_parstriangularis_sub-223-FUS','1';...
   'wleft_insula_sub-223-FUS','1';'wleft_lateralorbitofrontal_sub-223-FUS','1';'wleft_medialorbitofrontal_sub-223-FUS','1';'wleft_posteriorcingulate_sub-223-FUS','1';...
   'wleft_rostralanteriorcingulate_sub-223-FUS','1';'wleft_rostralmiddlefrontal_sub-223-FUS','1';'wleft_pericalcarine_sub-223-FUS', '1';'wright_Accumbens_sub-223-FUS','3';'wright_Amygdala_sub-223-FUS','2';...
   'wright_Caudate_sub-223-FUS','2';'wright_Pallidum_sub-223-FUS','2';'wright_Putamen_sub-223-FUS','2';'wright_Thalamus_sub-223-FUS','2';'wright_caudalanteriorcingulate_sub-223-FUS','1';...
   'wright_caudalmiddlefrontal_sub-223-FUS','1';'wright_frontalpole_sub-223-FUS','1';'wright_ifg_parsopercularis_sub-223-FUS','1';'wright_ifg_parsorbitalis_sub-223-FUS','1';'wright_ifg_parstriangularis_sub-223-FUS','1';...
   'wright_insula_sub-223-FUS','1';'wright_lateralorbitofrontal_sub-223-FUS','1';'wright_medialorbitofrontal_sub-223-FUS','1';'wright_posteriorcingulate_sub-223-FUS','1';'wright_rostralanteriorcingulate_sub-223-FUS','1';...
   'wright_rostralmiddlefrontal_sub-223-FUS','1'; 'wright_pericalcarine_sub-223-FUS','1'}%for cortical
counter = 1;
NAMES_CORTICAL = {};
NAME_ROICORTICAL = NAME_ROI(str2num(char(NAME_ROI(:,2))) ~= 2,:) 

for i = 1:size(subjects,1)
    for j = 1:3
                        if i < 10

        MAT = load([Path_results_Cortical 'resultsROI_Subject00' num2str(i) '_Condition00' num2str(j) '.mat']);%list of names for excel file
                        else
                                    MAT = load([Path_results_Cortical 'resultsROI_Subject0' num2str(i) '_Condition00' num2str(j) '.mat']);%list of names for excel file
                        end
        if i == 1 && j == 1

            for k = 1:size(MAT.names,2)
                NAMES_CORTICAL(counter) = NAME_ROICORTICAL(find(strcmp(MAT.names(k),NAME_ROICORTICAL(:,1))),1);

                counter = counter+1;
            end
        end
       NAMES_CORTICAL = erase( NAMES_CORTICAL,'_sub-223-FUS')

        INDEX_NAC_L = find(contains(MAT.names,'wleft_Accumbens_sub-223-FUS'));
        INDEX_NAC_R = find(contains(MAT.names,'wright_Accumbens_sub-223-FUS')); %atlas
        FC_Strength_Cortical(count,:) = MAT.Z(:,INDEX_NAC_R)';
        FC_Strength_Cortical(count+1,:) = MAT.Z(:,INDEX_NAC_L)';
        xlswrite(NAME,[{'ParticipantsName'},{'SESSION'}, {'SeedOfInterest'},NAMES_CORTICAL],'Cortical','A1:AE1')
        xlswrite(NAME,[repmat(string(subjects(i,:)),2,1) repmat(string(SESSION(j)),2,1) string({'Right_NAC';'Left_NAC'}) FC_Strength_Cortical(count:count+1,:) ],'Cortical',sprintf('A%s:AE%s',num2str(count+1),num2str(count+2)));
        count = count+2;

    end
end
count = 1;

%for subcortical
counter = 1;
NAMES_SUBCORTICAL = {};
NAME_ROISUBCORTICAL = NAME_ROI(str2num(char(NAME_ROI(:,2))) ~= 1,:)
for i = 1:size(subjects,1)
    for j = 1:3
                if i < 10

        MAT_SUBCORT = load([Path_results_Subcortical 'resultsROI_Subject00' num2str(i) '_Condition00' num2str(j) '.mat']);
                else
                            MAT_SUBCORT = load([Path_results_Subcortical 'resultsROI_Subject0' num2str(i) '_Condition00' num2str(j) '.mat']);
                end
       if i == 1 && j == 1
        for k = 1:size(MAT_SUBCORT.names,2)
            NAMES_SUBCORTICAL(counter) = NAME_ROISUBCORTICAL(find(strcmp(MAT_SUBCORT.names(k),NAME_ROISUBCORTICAL(:,1))),1);

            counter = counter+1;
        end
       end  
       NAMES_SUBCORTICAL = erase( NAME_ROISUBCORTICAL,'_sub-223-FUS')

       INDEX_NAC_L = find(contains(MAT_SUBCORT.names,'wleft_Accumbens_sub-223-FUS'));
        INDEX_NAC_R = find(contains(MAT_SUBCORT.names,'wright_Accumbens_sub-223-FUS')); %atlas
        FC_Strength_SubCortical(count,:) = MAT_SUBCORT.Z(:,INDEX_NAC_R)';
        FC_Strength_SubCortical(count+1,:) = MAT_SUBCORT.Z(:,INDEX_NAC_L)';
        xlswrite(NAME,[{'ParticipantsName'},{'SESSION'},{'SeedOfInterest'},NAMES_SUBCORTICAL(:,1)'],'Subcortical','A1:Z1')
        xlswrite(NAME,[repmat(string(subjects(i,:)),2,1) repmat(string(SESSION(j)),2,1) string({'Right_NAC';'Left_NAC'}) FC_Strength_SubCortical(count:count+1,:) ],'Subcortical',sprintf('A%s:Z%s',num2str(count+1),num2str(count+2)));
        count = count+2;

    end
end

end

% %% Posthoc bayesian with frontal area
% 
if FLAG_POSTHOC 
SESSION = {'sess-00','sess-07','sess-30'};
NAME ='Data_FS_NAC_POSTHOC.xlsx' %name of this excel file

count = 1;
%for frontal
counter = 1;
NAMES_POSTHOC = {};
for i = 1:size(subjects,1)
    for j = 1:3
        if i < 10
        MATFRONT = load([Path_Results_PostHoc 'resultsROI_Subject00' num2str(i) '_Condition00' num2str(j) '.mat']);%list of names for excel file
        else
                    MATFRONT = load([Path_Results_PostHoc 'resultsROI_Subject0' num2str(i) '_Condition00' num2str(j) '.mat']);%list of names for excel file
        end
                NAMES_POSTHOC = MATFRONT.names;

         

        INDEX_NAC_R = find(contains(MATFRONT.names,'wright_Accumbens_sub-223-FUS'));
        INDEX_NAC_L = find(strcmp(MATFRONT.names,'wleft_Accumbens_sub-223-FUS')); %atlas
        FC_Strength_Frontal(count,:) = MATFRONT.Z(:,INDEX_NAC_R)';
        FC_Strength_Frontal(count+1,:) = MATFRONT.Z(:,INDEX_NAC_L)';
        xlswrite(NAME,[{'ParticipantsName'},{'SESSION'}, {'SeedOfInterest'},NAMES_POSTHOC],'Frontal','A1:J1')
        xlswrite(NAME,[repmat(string(subjects(i,:)),2,1) repmat(string(SESSION(j)),2,1) string({'Right_NAC';'Left_NAC'}) FC_Strength_Frontal(count:count+1,:) ],'Frontal',sprintf('A%s:J%s',num2str(count+1),num2str(count+2)));
        count = count+2;
    end
    end

end
if FLAG_CONTROL
path_subjects='C:\Users\natha\Downloads\FUS_transferToNIDA(1)\FUS_transferToNIDA\';%To adapt accordingly
subjects=dir(path_subjects);
subjects={subjects.name}';
subjects = setdiff(subjects,{'.';'..'});
subjects = subjects(end-9:end,:);%To adapt accordingly to have the list of subjects to analyze

% brainRegions={'caudalanteriorcingulate','rostralanteriorcingulate', 'posteriorcingulate' , 'lateralorbitofrontal', ...
%    'medialorbitofrontal' ,'caudalmiddlefrontal'  , 'rostralmiddlefrontal' ,'frontalpole' ,'insula' ,'Thalamus','Caudate','Putamen','Pallidum','Amygdala','Accumbens'};



count = 1;
NAME_ROI = {'wleft_Accumbens_sub-223-FUS','3';'wleft_Amygdala_sub-223-FUS','2';'wleft_Caudate_sub-223-FUS','2';'wleft_Pallidum_sub-223-FUS','2';...
   'wleft_Putamen_sub-223-FUS','2';'wleft_Thalamus_sub-223-FUS','2';'wleft_caudalanteriorcingulate_sub-223-FUS','1';'wleft_caudalmiddlefrontal_sub-223-FUS','1';...
   'wleft_frontalpole_sub-223-FUS','1';'wleft_ifg_parsopercularis_sub-223-FUS','1';'wleft_ifg_parsorbitalis_sub-223-FUS','1';'wleft_ifg_parstriangularis_sub-223-FUS','1';...
   'wleft_insula_sub-223-FUS','1';'wleft_lateralorbitofrontal_sub-223-FUS','1';'wleft_medialorbitofrontal_sub-223-FUS','1';'wleft_posteriorcingulate_sub-223-FUS','1';...
   'wleft_rostralanteriorcingulate_sub-223-FUS','1';'wleft_rostralmiddlefrontal_sub-223-FUS','1'; 'wleft_pericalcarine_sub-223-FUS','3';'wright_Accumbens_sub-223-FUS','3';'wright_Amygdala_sub-223-FUS','2';...
   'wright_Caudate_sub-223-FUS','2';'wright_Pallidum_sub-223-FUS','2';'wright_Putamen_sub-223-FUS','2';'wright_Thalamus_sub-223-FUS','2';'wright_caudalanteriorcingulate_sub-223-FUS','1';...
   'wright_caudalmiddlefrontal_sub-223-FUS','1';'wright_frontalpole_sub-223-FUS','1';'wright_ifg_parsopercularis_sub-223-FUS','1';'wright_ifg_parsorbitalis_sub-223-FUS','1';'wright_ifg_parstriangularis_sub-223-FUS','1';...
   'wright_insula_sub-223-FUS','1';'wright_lateralorbitofrontal_sub-223-FUS','1';'wright_medialorbitofrontal_sub-223-FUS','1';'wright_posteriorcingulate_sub-223-FUS','1';'wright_rostralanteriorcingulate_sub-223-FUS','1';...
   'wright_rostralmiddlefrontal_sub-223-FUS','1'; 'wright_pericalcarine_sub-223-FUS','3'}%for cortical
counter = 1;
NAMES_CORTICAL = {};
NAME_ROICORTICAL = NAME_ROI(str2num(char(NAME_ROI(:,2))) ~= 2,:) 

for i = 1:size(subjects,1)
    for j = 1:3
         if i < 10
        MAT = load([Path_results_Cortical 'resultsROI_Subject00' num2str(i) '_Condition00' num2str(j) '.mat']);%list of names for excel file
        else
         MAT = load([Path_results_Cortical 'resultsROI_Subject0' num2str(i) '_Condition00' num2str(j) '.mat']);%list of names for excel file
         end
       
        if i == 1 && j == 1

            for k = 1:size(MAT.names,2)
                NAMES_CORTICAL(counter) = NAME_ROICORTICAL(find(strcmp(MAT.names(k),NAME_ROICORTICAL(:,1))),1);

                counter = counter+1;
            end
        end
NAME_ROICORTICAL = erase( NAME_ROICORTICAL(:,1),'_sub-223-FUS')
        INDEX_NAC_L= find(contains(MAT.names,'wleft_pericalcarine_sub-223-FUS'));
        INDEX_NAC_R = find(contains(MAT.names,'wright_pericalcarine_sub-223-FUS')); %atlas
        FC_Strength_Cortical(count,:) = MAT.Z(:,INDEX_NAC_R)';
        FC_Strength_Cortical(count+1,:) = MAT.Z(:,INDEX_NAC_L)';
        xlswrite(NAME,[{'ParticipantsName'},{'SESSION'}, {'SeedOfInterest'},NAMES_CORTICAL],'Control','A1:AE1')
        xlswrite(NAME,[repmat(string(subjects(i,:)),2,1) repmat(string(SESSION(j)),2,1) string({'LeftPericalcarine';'RightPericalcarine'}) FC_Strength_Cortical(count:count+1,:) ],'Control',sprintf('A%s:AE%s',num2str(count+1),num2str(count+2)));
        count = count+2;

    end
end
count = 1;

counter = 1;
NAMES_SUBCORTICAL = {};
NAME_ROISUBCORTICAL = NAME_ROI(str2num(char(NAME_ROI(:,2))) ~= 1,:) 

for i = 1:size(subjects,1)
    for j = 1:3
         if i < 10
        MAT = load([Path_results_Pericalcarine_Subcortical 'resultsROI_Subject00' num2str(i) '_Condition00' num2str(j) '.mat']);%list of names for excel file
        else
         MAT = load([Path_results_Pericalcarine_Subcortical 'resultsROI_Subject0' num2str(i) '_Condition00' num2str(j) '.mat']);%list of names for excel file
         end
       
        if i == 1 && j == 1

            for k = 1:size(MAT.names,2)
                NAMES_SUBCORTICAL(counter) = NAME_ROISUBCORTICAL(find(strcmp(MAT.names(k),NAME_ROISUBCORTICAL(:,1))),1);

                counter = counter+1;
            end
        end
NAME_ROISUBCORTICAL = erase( NAME_ROISUBCORTICAL(:,1),'_sub-223-FUS')
        INDEX_NAC_L= find(contains(MAT.names,'wleft_pericalcarine_sub-223-FUS'));
        INDEX_NAC_R = find(contains(MAT.names,'wright_pericalcarine_sub-223-FUS')); %atlas
        FC_Strength_SubCortical(count,:) = MAT.Z(:,INDEX_NAC_R)';
        FC_Strength_SubCortical(count+1,:) = MAT.Z(:,INDEX_NAC_L)';
        xlswrite(NAME,[{'ParticipantsName'},{'SESSION'}, {'SeedOfInterest'},NAMES_SUBCORTICAL],'ControlSubCortical','A1:AE1')
        xlswrite(NAME,[repmat(string(subjects(i,:)),2,1) repmat(string(SESSION(j)),2,1) string({'LeftPericalcarine';'RightPericalcarine'}) FC_Strength_SubCortical(count:count+1,:) ],'ControlSubCortical',sprintf('A%s:AE%s',num2str(count+1),num2str(count+2)));
        count = count+2;

    end
end
end