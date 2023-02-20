function [sig] = NormalizeSignal(sig, fs, varargin)
    % parse arguments  
    p = inputParser;
    p.addParameter('UsePlot',false);
    p.parse(varargin{:});
    
	UsePlot = p.Results.UsePlot;

    %% normalize the signal
    eupnea_secs = GetSection(sig, Pattern.EUPNEA);
    norm_sig = fillmissing(eupnea_secs(1).sig, 'linear');
    [ ~, ~, Y ] = ComputeSpectrum(norm_sig, fs);
    sig.sig = sig.sig./max(abs(Y));
    
    %% Plot result
    if UsePlot
        ref = ReferenceSig(fs);
        sigs = align_signals([sig, ref]);

        ShowSignal({sigs(:).sig}, fs);
        hold on
        plot([0 length(sig)/fs], [1 1], 'k');
        plot([0 length(sig)/fs], [-1 -1], 'k');
        hold off
    end
end

