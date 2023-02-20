function [comb_sig] = CombineSignals(sig1,sig2, ratio)
%COMBINESIGNALS Combines two signals with the given ratio

if nargin < 3
   ratio = 0.5; 
end

sig_len = min(length(sig1.sig), length(sig2.sig));

sig1.sig(sig1.labels == Pattern.UNDEFINED) = nan;
sig2.sig(sig2.labels == Pattern.UNDEFINED) = nan;
sig1.sig = fillmissing(sig1.sig, 'constant', 0); 
sig2.sig = fillmissing(sig2.sig, 'constant', 0); 

comb_sig.sig =  sig1.sig(1:sig_len,:).*ratio + ...
                sig2.sig(1:sig_len,:).*(1-ratio);
            
A_label = zeros(sig_len, 2);
A_label(:,1) = sig1.labels(1:sig_len);
A_label(:,2) = sig2.labels(1:sig_len);
A_label(A_label == Pattern.UNDEFINED) = nan;
comb_sig.labels = nanmedian(A_label, 2);
comb_sig.labels(isnan(comb_sig.labels)) = Pattern.UNDEFINED;

comb_sig.source = [num2str(ratio) '*' sig1.source '+' num2str(1-ratio) '*' sig2.source];
end

