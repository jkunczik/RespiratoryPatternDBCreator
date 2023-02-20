function [signal, info, ref_sig] = load_modalities(dir_name)
    folders = dir(dir_name);
    folders(~[folders.isdir]) = [];
    folders(1:2) = [];

    for i=1:length(folders)
        modality = upper(folders(i).name);
        [signal.(modality), info.(modality)] = load_signals(fullfile(dir_name, lower(modality)));
    end

    if isfield(signal, 'REF')
        ref_sig = signal.REF;
        signal = rmfield(signal, 'REF');
        info = rmfield(info, 'REF');
    end
end

