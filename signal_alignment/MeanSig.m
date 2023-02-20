function [sig] = MeanSig(sigs)
    % Computes the mean signal of all given signals

    % Crop signals to same length
    
    min_len = min(cellfun(@length, {sigs(:).sig}));
    A_sig = zeros(min_len, length(sigs));
    A_label = zeros(min_len, length(sigs));
    sig.source = 'mean of ';
    for i=1:length(sigs)
       A_sig(:,i) = sigs(i).sig(1:min_len);
       A_label(:,i) = sigs(i).labels(1:min_len);
       sig.source = [sig.source sigs(i).source ', '];

       A_sig(A_label(:,i) == Pattern.UNDEFINED, i) = nan; 
    end
   A_sig = fillmissing(A_sig, 'constant', 0);

    % Compute mean signal
    sig.sig = mean(A_sig, 2);

    A_label(A_label == Pattern.UNDEFINED) = nan;
    sig.labels = nanmedian(A_label, 2);
    sig.labels(isnan(sig.labels)) = Pattern.UNDEFINED;
end

