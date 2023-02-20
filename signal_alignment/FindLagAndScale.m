function [lag, p, q] = FindLagAndScale(orig_sig,ref_sig, fs, varargin)

    p = inputParser;
    p.addParameter('UsePlot',false);
    p.parse(varargin{:});

    UsePlot = p.Results.UsePlot;

    orig_sig = fillmissing(orig_sig, 'linear');
    ref_sig = fillmissing(ref_sig, 'linear');
    
    
    grid_len = 1000;
    max_corrs = nan(grid_len,1);
    max_lags = nan(grid_len, 1);
    f = linspace((fs-3),(fs+3),grid_len);
    
    for j=1:grid_len
        [p, q] = rat(f(j)/fs);           
        sig = resample(orig_sig, p, q);
        [corrs, lags] = xcorr(sig, ref_sig);
        [max_corrs(j), idx] = max(abs(corrs));
        max_lags(j) = lags(idx);
    end
    
    [~, idx] = max(max_corrs);
    lag = max_lags(idx);
    [p, q] = rat(f(idx)/fs);           
   
    if UsePlot
        sig = resample(orig_sig, p, q);
        ShowSignal({sig, [nan(max_lags(idx),1); ref_sig]}, fs);
    end
end

