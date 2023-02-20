function [sigs, info] = load_signals(data_dir)
    [~,root,~]=fileparts(data_dir);

    files = natsortfiles(dir(fullfile(data_dir, '*.mat')));
    for j = length(files):-1:1
        s = load(fullfile(files(j).folder, files(j).name));
        [~,name,~]=fileparts(files(j).name);

        if (size(s.sig, 2) > size(s.sig, 1))
            s.sig = s.sig';
        end
        s.labels = s.labels(:)';
        s.source = [root '/' name];
        
        sigs(j) = s;
    end

    infoStr = fileread(fullfile(data_dir, 'info.json'));
    info = jsondecode(infoStr);
end