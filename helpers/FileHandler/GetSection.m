function [sections] = GetSection(sig, label)

    sig.labels = sig.labels(:);
    if size(sig.sig, 2) > size(sig.sig, 1)
        sig.sig = sig.sig';
    end
    sig.labels(isnan(sig.labels)) = Pattern.UNDEFINED;
    index=find(sig.labels == label);
    if isempty(index)
        sections = {};
        return;
    end
    idx=find([index(1); diff(index)~=1; index(end)]); 

    if isempty(idx)
        lens = length(index);
    else
        lens=diff(idx);
    end   

    for i=length(lens):-1:1
        sections(i).source = sig.source;
        sections(i).start_id = index(idx(i));
    end

    sig_len = length(sig.sig);
    fn = fieldnames(sig);
    fn(contains(fn,{'labels', 'source'})) = [];
    for i=1:length(fn)
        if length(sig.(fn{i})) == sig_len
            C = mat2cell(sig.(fn{i})(index, :),lens, size(sig.(fn{i}), 2));
            for j=1:length(lens)
                sections(j).(fn{i}) = C{j};
            end
        end
    end

end

