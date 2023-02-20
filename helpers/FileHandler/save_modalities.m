function save_modalities(signal, info, dir)

    modalities = fieldnames(signal);

    for i=1:length(modalities)
        save_signals(fullfile(dir, lower(modalities{i})), signal.(modalities{i}), info);
    end
end

