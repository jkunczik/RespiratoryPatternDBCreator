function AlignDataSet(config)
%ALIGNDATASET Aligns all signals of the dataset with each other

[signal, info] = load_modalities(config.dirs.label_dir);

% 1 Align individual signals
% --------------------------
    modalities = fieldnames(signal);
    for i=1:length(modalities)
        sigs = signal.(modalities{i});
        inf = info.(modalities{i});

        for j=1:length(sigs)
            sigs(j) = NormalizeSignal(sigs(j), inf.fs);
            if contains(modalities{i}, {'BELT', 'RGB', 'CHEST_IRT'})
                sigs(j).sig = -sigs(j).sig;
            end
        end
        [sigs] = ResampleDataSet(sigs, inf, config.fs);
        sigs = align_signals(sigs);
        
        ShowSignal(cellfun(@(x) x(:,1), {sigs(:).sig}, 'UniformOutput', false), config.fs);

        signal.(modalities{i}) = sigs;
    end

% 2 Align signal classes with each other
% --------------------------------------
    med_sig_cell = cell(1, length(modalities) + 1);
    for i=1:length(modalities)
        med_sig.(modalities{i}) = NormalizeSignal(MeanSig(signal.(modalities{i})), config.fs);
        med_sig_cell{i} = med_sig.(modalities{i}).sig;
    end
    signal.REF = ReferenceSig(30);
    med_sig_cell{end} = signal.REF.sig;

    d = finddelays(med_sig_cell);

    for i=1:length(modalities)
        sigs = signal.(modalities{i});
        for j=1:length(sigs)
           sigs(j).sig = [nan(d(i), size(sigs(j).sig,2)); sigs(j).sig];
           sigs(j).labels = [nan(d(i), 1); sigs(j).labels];
           sigs(j).d = sigs(j).d + d(i);
        end
		signal.(modalities{i}) = sigs;
    end
    
%   Align Reference Signal
%   ~~~~~~~~~~~~~~~~~~~~~~
    signal.REF.sig = [nan(d(end), size(signal.REF.sig,2)); signal.REF.sig];
    signal.REF.labels = [nan(d(end), 1); signal.REF.labels];
    signal.REF.a_rr = [nan(d(end), 1); signal.REF.a_rr];
    signal.REF.rr = [nan(d(end), 1); signal.REF.rr];
    signal.REF.var_a_rr = [nan(d(end), 1); signal.REF.var_a_rr];
    signal.REF.var_rr = [nan(d(end), 1); signal.REF.var_rr];
    signal.REF.d = d(end);

%   Divide belt signal into chest and abdomen
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for i=length(signal.BELT):-1:1
       signal.CHEST_BELT(i) = signal.BELT(i);
       signal.CHEST_BELT(i).sig = signal.CHEST_BELT(i).sig(:,1);

       signal.ABD_BELT(i) = signal.BELT(i);
       signal.ABD_BELT(i).sig = signal.ABD_BELT(i).sig(:,2);
    end
    signal = rmfield(signal, 'BELT');
    modalities = fieldnames(signal);
    modalities(contains(modalities, 'REF')) = [];
    
%   verify result
%   ~~~~~~~~~~~~~
    names = cell(1, length(modalities)+1);
    for i=1:length(modalities)
        med_sig.(modalities{i}) = NormalizeSignal(MeanSig(signal.(modalities{i})), config.fs);
        med_sig_cell{i} = med_sig.(modalities{i}).sig;

        res = split(modalities{i}, '_');
        names{i} = '';
        for j=1:length(res)
            names{i} = [names{i} res{j}(1) lower(res{j}(2:end)) ' '];
        end
    end
    names{end} = 'Ref.';

    ShowSignal(med_sig_cell, config.fs);
    ShowSignal(signal.REF.sig, config.fs, 'NewFigure', false, 'PlotParams', {'Color', [0 0 0 0.5]});
    legend(names)

%   Save signals
%   ------------  
    save_modalities(signal, info, config.dirs.aligned_dir);
end

