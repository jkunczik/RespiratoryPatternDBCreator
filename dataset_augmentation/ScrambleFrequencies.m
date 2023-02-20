function [ds_sections] = ScrambleFrequencies(ds_sections, label, base_freq, freq_range, num_freqz)
    fn = fieldnames(ds_sections);
    fn(contains(fn, 'REF')) = [];

    sz = [length(ds_sections.(fn{1}).(label)), num_freqz];
    
    for i=sz(2):-1:1
        for j=sz(1):-1:1
           p = inf;
           q = inf;
           while (p*q > intmax('int32'))
               rr_new = round((diff(freq_range)*rand() + freq_range(1))*100)/100;
               rr_str = [' ' num2str(rr_new) ' BPM'];
               [p, q] = rat(ds_sections.CHEST_BELT.(label)(j).rr_med/rr_new);
           end
        
           idx = sub2ind(sz, j, i) + sz(1);
        
           for k=1:length(fn)
               ds_sections.(fn{k}).(label)(idx).sig = ...
                      resample(ds_sections.(fn{k}).(label)(j).sig, p, q);
               ds_sections.(fn{k}).(label)(idx).source = ...
                          [ds_sections.(fn{k}).(label)(j).source rr_str];
           end      
        end
    end
    % remove original signals
    for k=1:length(fn)
       ds_sections.(fn{k}).(label)(1:sz(1)) = [];
    end     
end