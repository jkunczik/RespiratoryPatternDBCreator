function [comb_sigs] = GenerateNewSignals(sigs, num_comb_sigs, fs)
%GENERATENEWSIGNALS Generate new signals by combining signals from the
%original dataset

    num_sigs = length(sigs);
    for i=num_comb_sigs:-1:1
                
        id_1 = randi(num_sigs);
        id_2 = randi(num_sigs);
        while (id_1 == id_2)
            id_2 = randi(num_sigs);
        end
        
        comb_sigs(i) = CombineSignals( ...
                                sigs(id_1), ...
                                sigs(id_2), ...
                                0.5*rand() + 0.25 ...
                            );
        comb_sigs(i) = NormalizeSignal(comb_sigs(i), fs);
        
    end
end

