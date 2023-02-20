function AugmentDataSet(config)
%AUGMENTDATASET Augments the given dataset by combining multiple signals
%together after randomly chaning their frequency

[signal, info, ref_sig] = load_modalities(config.dirs.aligned_dir);
    sig_names = fieldnames(signal);
    
    clear('sections');
    for k=1:length(sig_names)
        % generate new signals by combining original signals with random ratios
        comb_sig.(sig_names{k}) = GenerateNewSignals(signal.(sig_names{k}), config.num_comb_sigs, config.fs);

         % collect different patterns
        sections.(sig_names{k}) = GatherSections(comb_sig.(sig_names{k}));
        
        sec_names = fieldnames(sections.(sig_names{k}));
        % compute time characteristics
        for i=1:length(sec_names)
            sections = ComputeSectionFeatures(sections, sig_names{k}, sec_names{i}, config.fs);
        end
    end

    comb_sig.REF = ref_sig;
    sections.REF = GatherSections(ref_sig);
    
%   BRADYPNEA
%   ---------
%   10 random rates in range [5 10]
%   1 samples per recording
%   100 signal combinations
%   --> increase dataset size by a factor of 100
 sections = ScrambleFrequencies( ...
        sections, 'BRADYPNEA', ...
        5, [5 10], 10 ...
    );
  
%   EUPNEA
%   ------
%   5 random rates in range [12 18]
%   2 samples per recording
%   50 signal combinations
%   --> increase dataset size by a factor of 50
 sections = ScrambleFrequencies( ...
        sections, 'EUPNEA',...
        12, [12 18], 5 ...
    );
 sig_names = fieldnames(sections);
 for i=1:length(sig_names)
    for j=1:length(sections.(sig_names{i}).EUPNEA)
        norm_sig = ...
            fillmissing(sections.(sig_names{i}).EUPNEA(j).sig, 'linear');
        [ ~, ~, Y ] = ComputeSpectrum(norm_sig, config.fs);
        sections.(sig_names{i}).EUPNEA(j).sig = ...
            sections.(sig_names{i}).EUPNEA(j).sig./max(abs(Y));
    end
 end
  
%   HYPOPNEA
%   --------
%   10 random rates in range [12 18]
%   1 samples per recording
%   100 signal combinations
%   --> increase dataset size by a factor of 100
 sections = ScrambleFrequencies( ...
        sections, 'HYPOPNEA',...
        12, [12 18], 10 ...
    );

%   KUSSMAUL
%   -------------
%   10 random rates in range [20 35]
%   1 samples per recording
%   100 signal combinations
%   --> increase dataset size by a factor of 100
%   dismiss original Kussmaul sequency because of mismatch between demanded
%   and supplied effort by subjects. Use resampled hyperpnea instead
 modalities = fieldnames(sections);
 for i=1:length(modalities)
    sections.(modalities{i}).KUSSMAUL = sections.(modalities{i}).HYPERPNEA;
 end
 sections = ScrambleFrequencies( ...
        sections, 'KUSSMAUL',...
        10, [20 35], 10 ...
    );
  
%   HYPERPNEA
%   ---------
%   10 random rates in range [10 16]
%   1 samples per recording
%   10 signal combinations
%   --> increase dataset size by a factor of 100
 sections = ScrambleFrequencies( ...
        sections, 'HYPERPNEA',...
        10, [12 18], 10 ...
    );
  
%   TACHYPNEA
%   ---------
%   5 random rates in range [20 35]
%   2 samples per recording
%   100 signal combinations
%   --> increase dataset size by a factor of 100
 sections = ScrambleFrequencies( ...
        sections, 'TACHYPNEA', ...
        35, [20 35], 5 ...
    );

%   CHEYNE_STOKES
%   -------------
%   10 random rates in range [14 20]
%   1 samples per recording
%   100 signal combinations
%   --> increase dataset size by a factor of 100
 sections = ScrambleFrequencies( ...
        sections, 'CHEYNE_STOKES', ...
        20, [12 25], 10 ...
    );

 
%   BIOT
%   ----
%   5 random rates in range [14 20]
%   2 samples per recording
%   10 signal combinations
%   --> increase dataset size by a factor of 50
 sections = ScrambleFrequencies( ...
        sections, 'BIOT', ...
        15, [12 25], 5 ...
    );
  
%   APNEA
%   -----
%   10x Shorten or stretch the existing signals
%   10x Add noise with frequencies and amplitudes 
%   --> increase dataset size by a factor of 100
    fn = fieldnames(sections);
    fn(contains(fn, 'REF')) = [];
    sz = [length(sections.(fn{1}).APNEA), 9];
    
    for i=sz(2):-1:1
        for j=sz(1):-1:1

           l_new = (40*rand() + 20);
           len_str = [' ' num2str(l_new, 2) ' s'];

           t_new = (1:l_new*config.fs)/config.fs;

           idx = sub2ind(sz, j, i) + sz(1);
        
           for k=1:length(fn)
               t = 1:length(sections.(fn{k}).APNEA(j).sig);
               sections.(fn{k}).APNEA(idx).sig = ...
                      interp1(t, sections.(fn{k}).APNEA(j).sig, t_new)';
               sections.(fn{k}).APNEA(idx).source = ...
                          [sections.(fn{k}).APNEA(j).source len_str];
           end      
        end
    end
    mkdir(config.dirs.aug_ds_dir);
    save(fullfile(config.dirs.aug_ds_dir, 'sections.mat'), 'sections');


end

