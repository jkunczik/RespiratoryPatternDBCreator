
function [clust_sig] = ResampleDataSet(sigs, info, target_fs)

    if nargin < 3
        target_fs = info.fs;
    end 

    max_diff = 0.35;

    target_rr = 10;
    real_rrs = [];

    clust_sig = sigs;

    for i=1:length(sigs)
        section = GetSection(sigs(i), Pattern.HYPERPNEA);
        [~, ~, ~, ~, ~, rr_mean] = BreathFeatures(section.sig(:,1), info.fs, 'Prominence', 0.5, 'SpectrogramType', 'cwt');

        rr_mean(rr_mean == 0) = nan;
        real_rr = nanmedian(rr_mean);
        real_rrs(end+1) = real_rr;
        disp([sigs(i).source ': real_rr = ' num2str(real_rr)]);
    end

    % resample signals which don't have to much difference
    real_rr = median(real_rrs);
    [p, q] = rat(real_rr/target_rr*target_fs/info.fs);
    for i=1:length(clust_sig)
        sig = resample(clust_sig(i).sig, p, q);
        labels = ResampleLabels(clust_sig(i).labels, p, q)';
        sig_p = p;
        sig_q = q;

        if length(sig) > length(labels)
            sig(length(labels)+1:end,:) = [];
        elseif length(labels) > length(sig)
            labels(length(sig)+1:end) = [];
        end
        
        clust_sig(i).sig = sig;
        clust_sig(i).labels = labels;
        clust_sig(i).p = sig_p;
        clust_sig(i).q = sig_q;
    end

    al_sig = align_signals(clust_sig(abs(real_rr-real_rrs()) < max_diff));
    med_sig = NormalizeSignal(MeanSig(al_sig), target_fs);

    secs = GetSection(med_sig, Pattern.UNDEFINED);
    if (secs(1).start_id == 1)
        med_sig.sig(1:length(secs(1).sig)) = [];
        med_sig.labels(1:length(secs(1).sig)) = [];
    end

     for i=find(abs(real_rr-real_rrs) >= max_diff)
        disp(['high fs difference detected in sig: ' clust_sig(i).source '. resampling it individually']);

        [~, p, q] = FindLagAndScale(clust_sig(i).sig(:,1), med_sig.sig(:,1), target_fs);
        sig = resample(clust_sig(i).sig, p, q);
        labels = ResampleLabels(clust_sig(i).labels, p, q)';

        if length(sig) > length(labels)
            sig(length(labels)+1:end,:) = [];
        elseif length(labels) > length(sig)
            labels(length(sig)+1:end) = [];
        end

        clust_sig(i).sig = sig;
        clust_sig(i).labels = labels;
        clust_sig(i).p = p;
        clust_sig(i).q = q;
    end
end


