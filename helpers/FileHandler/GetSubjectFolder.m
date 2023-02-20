function [folders] = GetSubjectFolder(path)

    folders = struct2table(natsortfiles(dir(path)));
    folders = folders(folders.isdir,:);
    folders = folders(contains(folders.name, 'subj_'),:);
    folders = table2struct(folders);

end

