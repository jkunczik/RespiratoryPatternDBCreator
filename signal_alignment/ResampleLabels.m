function [new_labels] = ResampleLabels(labels,fs_new, fs_old)
    labels = labels(:);
    scale_fac = fs_new/fs_old;
    
    % find boundaries
    sig_bounds = [0; find(abs(diff(double(labels)))>=1); length(labels)];
    lens = sig_bounds(2:end) - sig_bounds(1:end-1);
    
    new_labels = [];
    for i=1:length(lens)
        label_sec = zeros(round(lens(i)*scale_fac), 1);
        label_sec(:) = labels(sig_bounds(i+1));
        new_labels = [new_labels; label_sec];
    end
end

