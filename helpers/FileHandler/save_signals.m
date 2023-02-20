function save_signals(data_dir, sigs, info)

    mkdir(data_dir)
   
    for j = 1:length(sigs)
        [~, name, ~] = fileparts(sigs(j).source);

        fn = fieldnames(sigs(j));
        fn(contains(fn,'source')) = [];

        for k=1:length(fn)
            eval([fn{k} '= sigs(j).' fn{k} ';']);
        end

        save(fullfile(data_dir, name),fn{:});
    end
    
    infoStr = jsonencode(info);
    
    fid = fopen(fullfile(data_dir, 'info.json'), 'w');
    fprintf(fid, infoStr);
    fclose(fid);
end