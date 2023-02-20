function ComputeSignalFeatures(data_dir, save_dir, config)
%COMPUTESIGNALFEATURES Summary of this function goes here
%   Detailed explanation goes here

    [signal, info, ref_sig] = load_modalities(data_dir);

%   Extract Features from whole signals
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sig_names = fieldnames(signal);
    for k=1:length(sig_names)
        for i=1:length(signal.(sig_names{k}))
             % remove invalid sections
             sig = signal.(sig_names{k})(i).sig;
             sig(signal.(sig_names{k})(i).labels == Pattern.UNDEFINED) = nan;
             sig = fillmissing(sig, 'constant', 0);

             [rr, a_rr, var_rr, var_a_rr, w_rr] = ...
                      BreathFeatures(sig, config.fs);
    
             signal.(sig_names{k})(i).w_rr = w_rr;
             signal.(sig_names{k})(i).rr = rr;
             signal.(sig_names{k})(i).a_rr = a_rr;
             signal.(sig_names{k})(i).var_rr = var_rr;
             signal.(sig_names{k})(i).var_a_rr = var_a_rr;
        end
    end

    signal.REF = ref_sig;

    save_modalities(signal, info, save_dir);
end

