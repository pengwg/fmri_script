%% Preprocessing FUS-OUD data
% This script is an adaptation of the preprocessing scripts made by
% Juliette Boscheron.
% I thank her for sharing that script with me.
% This script is preprocessing resting state fMRI data of patients (adults
% around 40 years old. It apply a voxel map displacement to correct for some
% potential distortion due to the magnetic field.
% To run this script you need: - structural file T1, functional file, a
% phase difference file and a magnitude file as well as the corresponding
% information (Phase direction encoding, echo time)
% IMPORTANT: PLease adapt the line indicated by a "%TO ADAPT ACCORDINGLY"

function Preprocessing_Function_Linux(subject_name, SESSIONS, TR, FWMH, blipdir, TotalReadOutTime, SPM_directives)

%% SETTINGS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Count_MotionFrame = 1;
%% SETTING UP VARIABLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% count = 1;

[path_spm, ~, ~] = fileparts(which('spm'));

for j = 1 : size(SESSIONS, 1)
    session_name = SESSIONS(j).name;

    path_to_func = [SESSIONS(j).folder '/' session_name '/func/'];
    path_to_anat = [SESSIONS(j).folder '/' session_name '/anat/'];
    path_to_fmap = [SESSIONS(j).folder '/' session_name '/fmap/'];
    path_to_ROI = [SESSIONS(j).folder '/' session_name '/ROI/'];

    %%  Get the functional scan
    FUNCFILE = [dir([path_to_func 'sub*task-RS*.nii']);
                dir([path_to_func '*task_rest*.nii']);
                dir([path_to_func '*_RS_-_FMRI*.nii'])];

    func_files = {FUNCFILE.name}';
    func_files = func_files(~contains(string(func_files),'P_A'),:); % so we remove from the list all the other type of functional file
    func_files = func_files(~contains(string(func_files),'PA'),:);
    func_files = func_files(~contains(string(func_files),'Cue'),:);
    func_files = strtrim(func_files);
    
    if size(func_files, 1) > 1
        disp('Warning: more than one candidate functional files found. Only process the first one.\n')
        disp(func_files)
    end

    func_file = func_files{1};

    % if size(func_files, 1)> 1 && contains(subject_name,'sub-222-FUS' )
    %     func_files = 'sub-222-FUS_ses-00_RS_-_FMRI_29.nii';
    % end
    % if size(func_files, 1)> 1 && contains(subject_name,'sub-220-FUS' )
    %     func_files = 'sub-220-FUS_ses-30_RS_-_FMRI_28.nii';
    % end

    % FILENAME(count) = {strtrim(FUNCFILE)} %to get the name of each file used
    % for each participant. Useful for debugging


    XX = spm_vol([path_to_func func_file]);
    %% Get number of slices
    Nb_slice = XX.dim; nb_slices = Nb_slice(3);
    nb_Frame = size(XX,1);
    %% Get voxel size
    %functional
    INFO = niftiinfo([path_to_func func_file]);
    INFO_Pix = INFO.PixelDimensions;
    voxel_size = [2.5 2.5 2.5];
    %structural
    voxel_size_T1 = [1 1 1];


    %% Get Slice timing value from json file
    json_func = regexprep(func_file, '\.nii$', '.json');
    json_data = jsondecode(fileread([path_to_func json_func]));
    slice_timing_ms = json_data.SliceTiming;

    %% Calculate VDM %%
    %%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_VDM
        ['****** Calculate VDM for subject ' subject_name '******* Session ' session_name  '********' ]
        clear matlabbatch

        fieldmap_phase = dir([path_to_fmap subject_name '*phasediff.nii']);
        json_fmap_mag = dir([path_to_fmap subject_name '*magnitude*.json']);
        
        EchoTime = [];
        for k = 1 : size(json_fmap_mag, 1)
            json_data = jsondecode(fileread([path_to_fmap json_fmap_mag(k).name]));
            EchoTime(k) = json_data.EchoTime;
        end

        [~, GOOD_FILE] = min(EchoTime); %give the minimum echotime from the 2 mag file (EchoTime_ToUse), and the index of where this file is relative to the list of json file
        fieldmap_mag = regexprep(json_fmap_mag(GOOD_FILE).name, '\.json$', '.nii');

        func_files_raw{j}=  [path_to_func func_file]; %spm_select('FPlist', path_to_func, '^sub.*/.nii$'); % this would work too
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {[path_to_fmap fieldmap_phase.name ',1']};
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {[path_to_fmap fieldmap_mag ',1']};
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [min(EchoTime)*1000 max(EchoTime)*1000];
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdir;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = TotalReadOutTime;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {[path_spm '/toolbox/FieldMap/T1.nii']};
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(1).epi = {[func_files_raw{j} ]};
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = 'g';
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
        spm_jobman('run', matlabbatch);
        clear matlabbatch

        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        func_files_vdm{j}=  spm_select('FPlist', path_to_fmap, '^vdm.*\.nii$');
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.scans = {[func_files_raw{j} ]};
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.data.vdmfile = {func_files_vdm{j}};
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.pedir = 2;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.which = [2 1];
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.rinterp = 4;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.mask = 1;
        matlabbatch{1}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'c';
        spm_jobman('run', matlabbatch);
    end
    clear fieldmap_mag fieldmap_phase func_files_raw matlabbatch

    %% Realign & Unwarp %%
    %%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_Resclice_Unwarp
        ['****** Realign & Unwarp for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        func_files_raw{j}=    [path_to_func  'c' func_file];

        %%
        func_files_vdm{j}=  spm_select('FPlist', [path_to_fmap], strcat(sprintf('^vdm.*session%d', j), '.*\.nii$'));

        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        matlabbatch{1}.spm.spatial.realignunwarp.data(1).scans = spm_file(cellstr(spm_select('expand', func_files_raw{j})));
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

        spm_jobman('run', matlabbatch);
    end
    clear func_files_vdm func_files_raw  matlabbatch

    %% Slice timing correction %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_SliceTiming_Correction

        ['****** Slice timing correction for subject' subject_name '******* Session' session_name '********' ]

        clear matlabbatch

        func_files_u = spm_select('FPlist', path_to_func, '^ucsub.*\.nii$'); % [path_to_func ls([path_to_func 'usub*task*.nii'])];%spm_select('FPlist', [path_to_func ], '^u.*/.nii$');
        matlabbatch{1}.spm.temporal.st.scans{1}= spm_file(cellstr(spm_select('expand', fullfile(func_files_u))));
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        matlabbatch{1}.spm.temporal.st.nslices = nb_slices;
        matlabbatch{1}.spm.temporal.st.tr = TR;
        matlabbatch{1}.spm.temporal.st.ta = TR-TR/nb_slices;
        matlabbatch{1}.spm.temporal.st.so = slice_timing_ms;
        matlabbatch{1}.spm.temporal.st.refslice = 0;
        matlabbatch{1}.spm.temporal.st.prefix = 'a';

        spm_jobman('run', matlabbatch);
    end

    clear func_files_u  matlabbatch slice_timing_ms


    %% Coregister anatomical to functional volumes %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_register_functional
        ['******  Coregister anatomical to functional volumes for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        mean_bold = spm_select('FPlist', path_to_func, '^meanc.*\.nii$');
        anat_vol_raw = spm_select('FPlist', path_to_anat, ['^' subject_name '.*\.nii$']);

        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {mean_bold};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[anat_vol_raw ',1']};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

        spm_jobman('run', matlabbatch);
    end
    clear mean_bold anat_vol_raw matlabbatch

    %% Segment anatomical scans in different tissue types %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_Segment
        ['****** Segment anatomical scans in different tissue types for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        anat_vol_registered = spm_select('FPlist', path_to_anat, ['^r' subject_name '.*\.nii$']);

        matlabbatch{1}.spm.spatial.preproc.channel.vols = {[anat_vol_registered ',1']};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[path_spm '/tpm/TPM.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[path_spm '/tpm/TPM.nii,2']};
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[path_spm '/tpm/TPM.nii,3']};
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[path_spm '/tpm/TPM.nii,4']};
        matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
        matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[path_spm '/tpm/TPM.nii,5']};
        matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
        matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[path_spm '/tpm/TPM.nii,6']};
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
        matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
        matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
            NaN NaN NaN];

        spm_jobman('run', matlabbatch);
    end
    clear anat_vol_registered matlabbatch

    %% Normalise functional volumes to standard space not slice corrected %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_normalize_functional
        ['******  Normalise functional volumes to standard space for subject ' subject_name '******* Session ' session_name '********NO SLICE TIMING CORRECTION*******' ]

        clear matlabbatch

        def_field = spm_select('FPlist', path_to_anat, '^y_rsub.*\.nii$');
        all_in_one_func_files_au = [];
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');
        func_files_au{j}=  spm_select('FPlist', [path_to_func ], '^uc.*\.nii$');

        volumes{j} = spm_file(cellstr(spm_select('expand', fullfile(func_files_au{j}))));
        all_in_one_func_files_au = cellstr([all_in_one_func_files_au; volumes{j}]);

        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {def_field};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = all_in_one_func_files_au;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = voxel_size;%TO ADAPT ACCORDINGLY
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        spm_jobman('run', matlabbatch);
    end
    clear def_field func_files_au volumes all_in_one_func_files_au matlabbatch

    %% Normalise functional volumes to standard space %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ['******  Normalise functional volumes to standard space for subject ' subject_name '******* Session ' session_name '******** SLICE TIMING CORRECTED ***********' ]

        clear matlabbatch
    if SPM_directives.do_normalize_functional

        def_field = spm_select('FPlist', path_to_anat, '^y_rsub.*\.nii$');
        all_in_one_func_files_au = [];
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        % if slice timing correction
        func_files_au{j}=  spm_select('FPlist', [path_to_func ], '^auc.*\.nii$');
        %if no slice timing correction
        % func_files_au{j}=  spm_select('FPlist', [path_to_func ], '^uc.*/.nii$');

        volumes{j} = spm_file(cellstr(spm_select('expand', fullfile(func_files_au{j}))));
        all_in_one_func_files_au = cellstr([all_in_one_func_files_au; volumes{j}]);

        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {def_field};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = all_in_one_func_files_au;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = voxel_size;%TO ADAPT ACCORDINGLY
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        spm_jobman('run', matlabbatch);
    end
    clear def_field func_files_au volumes all_in_one_func_files_au matlabbatch

    %% Normalise grey matter mask to standard space %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_register_greymatter
        ['******  Normalise grey matter mask to standard space for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        def_field = spm_select('FPlist', path_to_anat, '^y_rsub.*\.nii$');
        c1_rmask = spm_select('FPlist', path_to_anat, '^c1r.*\.nii$');
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');


        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {def_field};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {c1_rmask};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = voxel_size_T1;%TO ADAPT ACCORDINGLY
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        spm_jobman('run', matlabbatch);
    end

    clear def_field c1_rmask matlabbatch

    %% Normalise white matter mask to standard space %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_register_whitematter
        ['******  Normalise white matter mask to standard space for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        def_field = spm_select('FPlist', path_to_anat, '^y_rsub.*\.nii$');
        c2_rmask = spm_select('FPlist', path_to_anat, '^c2r.*\.nii$');
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');


        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {def_field};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {c2_rmask};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = voxel_size_T1;%TO ADAPT ACCORDINGLY
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        spm_jobman('run', matlabbatch);
    end
    clear def_field c2_rmask matlabbatch

    %% Normalise csf mask to standard space %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_register_csf
        ['******  Normalise csf mask to standard space for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        def_field = spm_select('FPlist', path_to_anat, '^y_rsub.*\.nii$');
        c3_rmask = spm_select('FPlist', path_to_anat, '^c3r.*\.nii$');
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');


        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {def_field};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {c3_rmask};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = voxel_size_T1;%TO ADAPT ACCORDINGLY
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        spm_jobman('run', matlabbatch);
    end
    clear def_field c3_rmask matlabbatch

    %% Normalise structural T1 to standard space %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_register_structural
        ['******  Normalise structural T1 to standard space for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        def_field = spm_select('FPlist', path_to_anat, '^y_rsub.*\.nii$');
        struct = spm_select('FPlist', path_to_anat, '^sub.*\.nii$');
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');


        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {def_field};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {struct};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = voxel_size_T1;%TO ADAPT ACCORDINGLY
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

        spm_jobman('run', matlabbatch);
    end
    clear def_field struct matlabbatch

    %% Smooth functional volumes %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_smooth

        ['******  Smooth functional volumes for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        all_in_one_func_files_wau = [];
        % if slice timing correction
        %  func_files_wau=  spm_select('FPlist', [path_to_func ], '^wau.*/.nii$');
        % if no slice timing correction
        func_files_wau=  spm_select('FPlist', path_to_func, '^wu.*\.nii$');

        volumes = spm_file(cellstr(func_files_wau));
        all_in_one_func_files_wau = cellstr([all_in_one_func_files_wau; volumes]);
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');


        matlabbatch{1}.spm.spatial.smooth.data = all_in_one_func_files_wau;
        matlabbatch{1}.spm.spatial.smooth.fwhm = [FWMH FWMH FWMH];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';

        spm_jobman('run', matlabbatch);
    end
%     Flag_smooth = 1;
    clear all_in_one_func_files_wau matlabbatch

        %% Smooth functional volumes 8 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SPM_directives.do_smooth

        ['******  Smooth functional volumes for subject ' subject_name '******* Session ' session_name '********' ]

        clear matlabbatch

        all_in_one_func_files_wau = [];
        % if slice timing correction
        %  func_files_wau=  spm_select('FPlist', [path_to_func ], '^wau.*/.nii$');
        % if no slice timing correction
        func_files_wau=  spm_select('FPlist', [path_to_func ], '^wu.*\.nii$');

        volumes = spm_file(cellstr(func_files_wau));
        all_in_one_func_files_wau = cellstr([all_in_one_func_files_wau; volumes]);
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');


        matlabbatch{1}.spm.spatial.smooth.data = all_in_one_func_files_wau;
        matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's8';

        spm_jobman('run', matlabbatch);
    end
%     Flag_smooth = 1;
    clear all_in_one_func_files_wau matlabbatch

    %% Smooth functional volumes %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ['******  Smooth functional volumes for subject ' subject_name '******* Session ' session_name '******** with slice timing correction' ]

    clear matlabbatch

    all_in_one_func_files_wau = [];
    func_files_wau=  spm_select('FPlist', [path_to_func ], '^wauc.*\.nii$');
    volumes = spm_file(cellstr(func_files_wau));
    all_in_one_func_files_wau = cellstr([all_in_one_func_files_wau; volumes]);
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');


    matlabbatch{1}.spm.spatial.smooth.data = all_in_one_func_files_wau;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [FWMH FWMH FWMH];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

    spm_jobman('run', matlabbatch);

    clear all_in_one_func_files_wau matlabbatch

    %% END %%

    %% SANITY CHECK: FRAMEWISE DISPLACEMENT
    %Get the number of frame that have a displacement bigger than a certain
    %threshold (standard is the size of a voxel and save a picture of the displacement in mm in the three axis
    if SPM_directives.do_MotionCheck
        [Number_Of_Frame_Above_Threshold_Per_Session] = Check_motionDisplacement_NHM_v1(path_to_func, 0.5, session_name, subject_name);
        PercentDisplacement(Count_MotionFrame) = Number_Of_Frame_Above_Threshold_Per_Session*100/nb_Frame;
        Subject_Session(Count_MotionFrame) = {[subject_name '-' session_name]}
        Count_MotionFrame = Count_MotionFrame + 1;
    end
end
%reset each flag to zero after each session

end
