
% It is expected that the eval_dir contains at least one folder
% "1_extracted_signals", containing the extracted raw signals. This folder
% can be downloaded from PhysioNet. 
eval_dir = '.';

config.fs = 30;
config.num_comb_sigs = 100;

config.sigs(1).name = 'border_irt';
config.sigs(1).file_name = 'border_irt.mat';
config.sigs(1).fs = 30;

config.sigs(2).name = 'chest_irt';
config.sigs(2).file_name = 'chest_irt.mat';
config.sigs(2).fs = 30;

config.sigs(3).name = 'rgb';
config.sigs(3).file_name = 'rgb_fused.mat';
config.sigs(3).fs = 15;

config.sigs(4).name = 'belt';
config.sigs(4).file_name = 'belt.mat';
config.sigs(4).fs = 30;

config.dirs.extracted_sigs_folder = fullfile(eval_dir, '1_extracted_signals');
config.dirs.label_dir = fullfile(eval_dir, '2_labeled_signals');
config.dirs.aligned_dir = fullfile(eval_dir, '3_aligned_signals');
config.dirs.orig_ds_dir = fullfile(eval_dir, '4_original_dataset');
config.dirs.aug_ds_dir = fullfile(eval_dir, '5_augmented_dataset');
config.dirs.wfdb_dir = fullfile(eval_dir, '6_wfdb');


start_from_step = 5;

%% (1) Label and normalize Data
%  ============================

if (start_from_step <= 1)
    LabelDataSet(config)
end

%% (2) Align Signals
%  ================

if (start_from_step <= 2)
    AlignDataSet(config)
end

%% (3) Compute features from original dataset
%  ==========================================
if (start_from_step <= 3)
    ComputeSignalFeatures(config.dirs.aligned_dir, config.dirs.orig_ds_dir, config);
    
    % Gather sections from full length signals
    [signal, info, ref_sig] = load_modalities(config.dirs.orig_ds_dir);
    sig_names = fieldnames(signal);
    for k=1:length(sig_names)
        sections.(sig_names{k}) = GatherSections(signal.(sig_names{k}));
    end
    sections.REF = GatherSections(ref_sig);
    save(fullfile(config.dirs.orig_ds_dir, 'sections.mat'), 'sections');

    ComputeAllSectionFeatures(config.dirs.orig_ds_dir, config);
end

%% (4) Augment Dataset
%  ===================
if (start_from_step <= 4)
    AugmentDataSet(config)
end

%% (5) Extract Features from augmented dataset
%  ===========================================
if (start_from_step <= 5)
    ComputeAllSectionFeatures(config.dirs.aug_ds_dir, config)
end

%% (7) Export to WFDB format
%  =========================
if (start_from_step <= 7)
    curr_dir = pwd;
    cd(wfdb_dir);

    [signal, info, ref_sig] = load_modalities(orig_ds_dir);
    modalities = fieldnames(signal);
    % ATTENTION: This file and folder are not provided with the dataset,
    % due to privacy concerns.
    subj_info = readtable(fullfile(video_folder, "Probanden.xlsx"));

    % Raw signals
    subjects = natsortfiles(dir(fullfile(extracted_sigs_folder, 'subj_*')));
    clear('raw_signal');
    for i=1:length(modalities)
        raw_signal.(modalities{i}).sig = [];
    end
    for i=1:length(subjects)
        belt_file = fullfile(subjects(i).folder, subjects(i).name, 'belt.mat');
        if isfile(belt_file)
            load(belt_file);
            raw_signal.CHEST_BELT(end+1).sig = res(:,1);
            raw_signal.ABD_BELT(end+1).sig = res(:,2);
        end
        rgb_file = fullfile(subjects(i).folder, subjects(i).name, 'rgb_raw.mat');
        if isfile(rgb_file)
            load(rgb_file);
            raw_signal.RGB(end+1).sig = GatherSignalComponents(posi, @SVDMax);
        end
        chest_irt_file = fullfile(subjects(i).folder, subjects(i).name, 'chest_irt.mat');
        if isfile(chest_irt_file)
            load(chest_irt_file);
            raw_signal.CHEST_IRT(end+1).sig = GatherSignalComponents(posi, @SVDMax);
        end
        border_irt_file = fullfile(subjects(i).folder, subjects(i).name, 'border_irt.mat');
        if isfile(border_irt_file)
            load(border_irt_file);
            raw_signal.BORDER_IRT(end+1).sig = sig;
        end
    end
    for i=1:length(modalities)
        raw_signal.(modalities{i})(1) = [];
    end

    for i=1:length(modalities)
        modality = modalities{i};
        save_dir = fullfile('raw', lower(modality));
        mkdir(save_dir);
        pushd(save_dir);

        for j=1:length(raw_signal.(modality))
            sig = signal.(modality)(j);
            raw_sig = raw_signal.(modality)(j);
            raw_sig.sig = ResamplePQ(raw_sig.sig, sig.p, sig.q);
            raw_sig.sig = [nan(sig.d, size(raw_sig.sig, 2)); raw_sig.sig];

            [~, name, ~] = fileparts(signal.(modality)(j).source);
            res = split(name, '_');
            number_str = res{2};
            number = str2double(number_str);

            X = raw_sig.sig;
            X(:, all(isnan(raw_sig.sig))) = [];
            fname = number_str;
            bit_res = 16;
            adu = [];
            info = {['<source>: ' name],...
                    ['<age>: ' num2str(subj_info{number,'Alter'}) '  <sex>: ' subj_info{number,'Geschlecht'}{1} '  <diagnoses>: (none)  <medications>: (none)']};
            gain = [];
            sg_name = [];
        
            [xbit] = mat2wfdb(X,fname,lwir_fs,bit_res,adu,info,gain,sg_name);
        end
        popd
    end

% Reference signal
    save_dir = fullfile('reference');
    mkdir(save_dir);
    pushd(save_dir);

    sig = ref_sig;
    X = [sig.sig, sig.rr, sig.a_rr, sig.var_rr, sig.var_a_rr];
    fname = 'reference';
    bit_res = 16;
    adu = {'n.u', 'Breaths/min', 'n.u.', 'Breaths/min', 'n.u.'};
    info = {'<source>: reference'};
    gain = [];
    sg_name = {'Resp', 'RR', 'A_RR', 'RR_var', 'A_RR_var'};
    
    [xbit] = mat2wfdb(X,fname,lwir_fs,bit_res,adu,info,gain,sg_name);
    popd

    % Extracted signals
    load(fullfile(orig_ds_dir, 'sections.mat'));

    for i=1:length(modalities)
        modality = modalities{i};
        save_dir = fullfile('extracted', lower(modality));
        mkdir(save_dir);
        pushd(save_dir);

        % signals 
        for j=1:length(signal.(modality))
            sig = signal.(modality)(j);
            [~, name, ~] = fileparts(signal.(modality)(j).source);
            res = split(name, '_');
            number_str = res{2};
            number = str2double(number_str);

            X = [sig.sig, sig.labels', sig.rr, sig.a_rr, sig.w_rr, sig.var_rr, sig.var_a_rr];
            fname = number_str;
            bit_res = 16;
            adu = {'n.u', '', 'Breaths/min', 'n.u.', 's', 'Breaths/min', 'n.u.'};
            info = {['<source>: ' name],...
                    ['<age>: ' num2str(subj_info{number,'Alter'}) '  <sex>: ' subj_info{number,'Geschlecht'}{1} '  <diagnoses>: (none)  <medications>: (none)']};
            gain = [];
            sg_name = {'Resp', 'label', 'RR', 'A_RR', 'W_RR', 'RR_var', 'A_RR_var'};
        
            [xbit] = mat2wfdb(X,fname,lwir_fs,bit_res,adu,info,gain,sg_name);
        end
        popd
        
        %sections
        sec_names = fieldnames(sections.(modality));
        for k=1:length(sec_names)
            section = sections.(modality).(sec_names{k});
            sec_dir = fullfile(save_dir, lower(sec_names{k}));
            mkdir(sec_dir);
            pushd(sec_dir);

            for j=1:length(section)
                sig = section(j);
                X = [sig.sig, sig.rr, sig.a_rr, sig.w_rr, sig.var_rr, sig.var_a_rr];
                fname = num2str(j);
                bit_res = 16;
                adu = {'n.u', 'Breaths/min', 'n.u.', 's', 'Breaths/min', 'n.u.'};
                info = {['<source>: ' sig.source]};
                gain = [];
                sg_name = {'Resp', 'RR', 'A_RR', 'W_RR', 'RR_var', 'A_RR_var'};
            
                [xbit] = mat2wfdb(X,fname,lwir_fs,bit_res,adu,info,gain,sg_name);
            end
            stats = table([section.rr_med]', [section.rr_var]', [section.a_med]', [section.a_var]', 'VariableNames', {'RR_med', 'RR_var', 'A_med', 'A_var'});
            writetable(stats,'stats.csv')
            popd;
        end
    end

    % augmented
    load(fullfile(aug_ds_dir, 'sections.mat'));
    modalities = fieldnames(sections);
    modalities(contains(modalities, 'REF')) = [];

    for i=1:length(modalities)
        modality = modalities{i};
        save_dir = fullfile('augmented', lower(modality));
        mkdir(save_dir);
        
        %sections
        sec_names = fieldnames(sections.(modality));
        for k=1:length(sec_names)
            section = sections.(modality).(sec_names{k});
            sec_dir = fullfile(save_dir, lower(sec_names{k}));
            mkdir(sec_dir);
            pushd(sec_dir);

            for j=1:length(section)
                sig = section(j);
                X = [sig.sig, double(sig.rr)/1e3, double(sig.a_rr)/1e3, double(sig.w_rr)/1e3, double(sig.var_rr)/1e2, double(sig.var_a_rr)/1e5];
                fname = num2str(j);
                bit_res = 16;
                adu = {'n.u', 'Breaths/min', 'n.u.', 's', 'Breaths/min', 'n.u.'};
                info = {['<source>: ' sig.source]};
                gain = [];
                sg_name = {'Resp', 'RR', 'A_RR', 'W_RR', 'RR_var', 'A_RR_var'};
            
                [xbit] = mat2wfdb(X,fname,lwir_fs,bit_res,adu,info,gain,sg_name);
            end
            sec = sections.(modality).(sec_names{k});
            stats = table([sec.rr_med]', [sec.rr_var]', [sec.a_med]', [sec.a_var]', 'VariableNames', {'RR_med', 'RR_var', 'A_med', 'A_var'});
            writetable(stats,'stats.csv')
            popd;   
        end
    end

    % create RECORDS
    subdirs = create_wfdb_records(pwd);

    root_path = regexprep(pwd,'\\','\\\\');

    fid = fopen(fullfile(root_path, 'RECORDS'), "w");
    for i=1:length(subdirs)
        sub_path = regexprep(subdirs{i}, root_path, '');
        sub_path = regexprep(sub_path, '\\', '/');
        sub_path = regexprep(sub_path, '^/', '');

        fprintf(fid, '%s/\n', sub_path);
    end
    fclose(fid);

    cd(curr_dir);
end

%% Sub Functions
%  =============

function [] = plotmatrix_section(sections, modality)
    sec_names = fieldnames(sections.(modality));
    sec_len = length(sections.(modality).(sec_names{1}));
    X = nan(sec_len*length(sec_names), 4);
    Label = cell(sec_len*length(sec_names), 1);
    for i=1:length(sec_names)
        sec = sections.(modality).(sec_names{i});
        X((1:sec_len)+(i-1)*sec_len, :) = ...
            [[sec.rr_med]', [sec.a_med]', [sec.rr_var]', [sec.a_var]'];
        Label((1:sec_len)+(i-1)*sec_len) = {sec_names{i}};
    end
    colors = lines(length(sec_names));
    xnames = {'rr_{med}', 'a_{med}', 'rr_{var}', 'a_{var}'};
    gplotmatrix(X,[], Label, colors, [], 5,[],'grpbars',xnames);
end

function dirs = create_wfdb_records(path)
    files = natsortfiles(dir(fullfile(path, '*.hea')));
    dirs = {};
    if ~isempty(files)
        fid = fopen(fullfile(path, 'RECORDS'), "w");
        for i=1:length(files)
            [~, filename] = fileparts(files(i).name);
    
            fprintf(fid, '%s\n', filename);
        end
        fclose(fid);
        dirs{end+1} = path;
    end

    sub_dirs = dir(path);
    sub_dirs(~[sub_dirs.isdir]) = []; 
    sub_dirs(matches({sub_dirs.name} ,{'.','..'})) = [];
    for i = 1 : length(sub_dirs)
      subdir = sub_dirs(i).name;
      dirs = [dirs; create_wfdb_records(fullfile(path, subdir))];
    end
end