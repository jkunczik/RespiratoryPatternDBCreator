function [sigs, err] = align_signals(sigs)
    % remove invalid sections with nans
    cleaned_sigs = cell(length(sigs), 1);
    for i = 1:length(sigs)
        if size(sigs(i).sig, 2) > size(sigs(i).sig,1)
            sigs(i).sig = sigs(i).sig';
        end
        sigs(i).labels = sigs(i).labels(:);

        cleaned_sigs{i} =  sigs(i).sig;
        cleaned_sigs{i}(sigs(i).labels == Pattern.UNDEFINED) = nan; 
        cleaned_sigs{i} = fillmissing(sigs(i).sig, 'constant', 0); 
    end
    % align signals
    [d, err] = finddelays(cleaned_sigs);
    
    for i=1:length(sigs)
        sigs(i).sig = [nan(d(i), size(sigs(i).sig,2)); sigs(i).sig];
        sigs(i).labels = [nan(d(i), 1); sigs(i).labels];
        sigs(i).d = d(i);
    end
end