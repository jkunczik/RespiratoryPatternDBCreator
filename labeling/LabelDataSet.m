function [] = LabelDataSet(config)
%LABELDATASET Labels all raw signals of the data set.

    folders = GetSubjectFolder(config.dirs.extracted_sigs_folder);

    for i=1:length(folders)
        for j=1:length(config.sigs)
            try
                info = config.sigs(j);
                disp(['Labeling signal ' info.name ' from folder ' folders(i).name])

                content = load(fullfile(folders(i).folder, folders(i).name, info.file_name));    
                labels = LabelSignal(content.res, info.fs);
                sig.sig = detrend(content.res);
                sig.labels = labels;
                sig.source = folders(i).name;
        
                sig = NormalizeSignal(sig, info.fs);
                save_signals(fullfile(config.dirs.label_dir, info.name), sig, info);
            catch e
                warning(['Cloud not label signal signal' info.name ' in folder ' folders(i).name ' . Reason: ' e.message]);
            end
        end
    end
end

