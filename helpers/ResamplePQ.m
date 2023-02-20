function [fixedSig] = ResamplePQ(y, p, q)

pad = 100;

y_pre = flipud(y(1:pad,:));
y_post = flipud(y(end-pad+1:end,:));

fixedSig = resample([y_pre; y; y_post], p, q);
padpoints = round(pad*p/q);
fixedSig = fixedSig(padpoints:end-padpoints, :);

end

