function [rr, a_rr, var_rr, var_a_rr, w_rr, rr_mean] = BreathFeatures(sig,fs, varargin)

p = inputParser;
p.addParameter('upfac',5);
p.addParameter('Prominence', 0.3);
p.addParameter('UsePlot',false);
p.addParameter('AutoClose',true);
p.addParameter('RR_max', 40);
p.addParameter('RR_min', 3);
p.addParameter('Normalize', false);
p.addParameter('SpectrogramType', 'cwt'); %'fft' 'wsst' or 'cwt'
p.parse(varargin{:});

upfac = p.Results.upfac;
Prominence = p.Results.Prominence;
UsePlot = p.Results.UsePlot;
AutoClose = p.Results.AutoClose;
RR_max = p.Results.RR_max;
RR_min = p.Results.RR_min;
Normalize = p.Results.Normalize;
SpectrogramType = p.Results.SpectrogramType;

window = 'hamming';
win_length = 30;
spec_fs = 3;

min_amplitude = 0.15;
rr_range = [4 40];
rr_eps = 0.05;
max_rr_delta = 5;
jump_eps = 9; %s

findpeaks_config = {
        'MinPeakProminence',Prominence, ...
        'WidthReference','halfheight', ...
        'MaxPeakWidth', 2/3*60/RR_min, ...
        'MinPeakDist', 60/RR_max,...
    };

spectrogram_config = {
    'WinLength', win_length, ...
    'Window', window, ...
    'MinFreqRes', 0.005, ...
    'Type', 'imagesc', ...
    'Scale', 'logarithmic',...
    'SpectrogramType', SpectrogramType,...
};

sig = sig(:);
sig_len = length(sig);


%% preprocess signal
% remove nan and signal boundaries
ids = find(~isnan(sig));
sig = sig(ids(1) : ids(end));
t_start = (min(ids)-1)/fs; % in s

% fill missing sections with 0
sig = fillmissing(sig, 'constant', 0);

% smooth signal
sig = sgolayfilt(sig, 2, round(1.5*fs));

% remove spikes
% spike_ids = find(abs(sig) > 4);
% large_sig_mask = false(length(sig), 1);
% for i=2:length(spike_ids)
%     if (spike_ids(i)/fs < spike_ids(i-1)/fs + jump_eps)
%         large_sig_mask(spike_ids(i-1):spike_ids(i)) = 1;
%     end
% end
% ids = [find(abs(diff(large_sig_mask))>0); length(sig)];
% large_sig = ~large_sig_mask(ids);
% for i=1:length(ids)-1
%     range = [ids(i), (ids(i+1))];
%     if large_sig(i)
%         sec_len = diff(range)/fs;
%         if  sec_len < jump_eps
%             sig(range(1):range(2)) = nan;
%         end
%     end
% end


%% pad signal
% % 1. find first two peaks
% [~, u_locs] = findpeaks(sig, fs, 'NPeaks', 3, findpeaks_config{:});
% [~, l_locs] = findpeaks(-sig, fs, 'NPeaks', 3, findpeaks_config{:});
% 
% locs = sort([u_locs; l_locs]);
% 
% % 2. determine necessary pad_length
% if length(locs) >= 4 
%     pad_length = ceil(max([diff(locs([2,4])), win_length/2])); % in s
% else
%     pad_length = win_length/2;
% end
% 
% if length(locs) < 1
%     locs(1) = 1/fs;
% end
% 
% % 3. pad signal by mirroring with pad_length
% if length(locs) < 1 ||  (2*locs(1) + pad_length)*fs > length(sig)
%     sig = [ sig(1)*ones(pad_length*fs,1);
%             sig
%             sig(end)*ones(pad_length*fs,1)];
% else  
%     sig = [flip(sig(round(locs(1)*fs + (1:(pad_length+locs(1))*fs))));
%            sig(locs(1)*fs:end)
%            sig(end)*ones(pad_length*fs,1)];
% end
pad_length = 1/fs;

%% normalize the signal
if Normalize
    ShowSignal(sig, fs, 'FullScreen', true);
    [norm_range, ~] = ginput(2);
    ids = round(fs*[norm_range(1) : norm_range(2)]);
    close(gcf);

    norm_sig = sig(ids);
    [ ~, f, Y ] = ComputeSpectrum(norm_sig, fs);
    sig = sig*1/max(2*abs(Y));
end

%% Compute temporal characteristics
up_sig = ResamplePQ(sig, upfac, 1);

% subtract lower envelope
% [peaks, locs] = findpeaks(-up_sig, 'MinPeakProminence',Prominence);
% peaks = [up_sig(1); peaks; up_sig(end)];
% locs = [1; locs; length(up_sig)]/(fs*upfac);
% env = interp1(locs, -peaks, locs(1):1/(fs*upfac):locs(end), 'linear');
% 
% up_sig = up_sig - env';

[~, locs, w, p] = findpeaks(up_sig, fs*upfac, findpeaks_config{:});
locs = [0; locs; length(sig)/fs];
p = [0; p; 0];
w = [0; w; 0];

t = locs;
w_rr_t = w;
a_rr_inst_t = p/2;
rr_inst_t = [0; 1./diff(locs)*60];

% Filter values with to low or high frequency
rr_inst_t(rr_inst_t < rr_range(1) + rr_eps) = nan;
rr_inst_t(rr_inst_t > rr_range(2) - rr_eps) = nan;
a_rr_inst_t(isnan(rr_inst_t)) = 0;

%% Compute spectral characteristics
spec_sig = ResamplePQ(fillmissing(sig, 'const', 0), spec_fs, fs);
if strcmpi(SpectrogramType, 'fft')
    [P, f, t] = ComputeSpectrogram(spec_sig, spec_fs, spectrogram_config{:});
    snr_down = nan(size(P,2), 1);
    for i = 1:size(P,2)
        snr_down(i) = snr(P(:,i), f, 3, 'psd');
    end
else
    if strcmpi(SpectrogramType, 'cwt')
        [P, f] = cwt(spec_sig, spec_fs, 'TimeBandwith', 10, 'VoicesPerOctave', 48, 'FrequencyRange', rr_range/60);
    elseif strcmpi(SpectrogramType, 'wvd')
        [P,f,t] = wvd(spec_sig, spec_fs, 'smoothedPseudo');
    else
        [P, f] = wsst(spec_sig, spec_fs, 'VoicesPerOctave', 20);
        P(f < rr_range(1)/60 | f > rr_range(2)/60,:) = [];
        f(f < rr_range(1)/60 | f > rr_range(2)/60) = [];
        f = f(:);
    end
    P = abs(P);
    t_mean = (0:(size(P,2)-1))'/spec_fs;
    snr_down = 10 * log10(max(P)/0.3);
end

[a_rr_mean_t, idx] = max(P);
a_rr_mean_t = a_rr_mean_t';
rr_mean_t = f(idx)*60;


%% Artifact correction and resampling
t_res = (t(1):1/fs:t(end))';
start_pad = nan(round(t(1)*fs),1);
method = 'next';

% Temporal
a_rr_mean_res = [start_pad; interp1(t_mean, a_rr_mean_t, t, method)];
if length(a_rr_mean_res) >= length(rr_inst_t)
    a_rr_mean_res(length(rr_inst_t)+1:end) = [];
else
    a_rr_mean_res = [a_rr_mean_res;
                a_rr_mean_res(end)*ones((length(rr_mean_t)-length(a_rr_mean_res)), 1)];
end
low_amp_mask =  a_rr_inst_t <= 0.05 | a_rr_mean_res  <= 0.05;
rr_inst_t = ArtifactCorrection(rr_inst_t, low_amp_mask, t);

a_rr_inst_t(rr_inst_t == 0) = 0;

rr_inst = [start_pad; interp1(t, rr_inst_t, t_res, method)];
% rr_inst = rr_inst(1:end-1);

a_rr_inst = [start_pad; interp1(t, a_rr_inst_t, t_res, method)];
w_rr = [start_pad; interp1(t, w_rr_t, t_res, method)];

tail = ones(length(sig) - length(rr_inst),1);

rr_inst = [rr_inst; tail*rr_inst(end)];
a_rr_inst = [a_rr_inst; tail*a_rr_inst(end)];
w_rr = [w_rr; tail*w_rr(end)];

% Spectral
a_rr_res = [start_pad; interp1(t, a_rr_inst_t, t_mean, method)];
if length(a_rr_res) >= length(rr_mean_t)
    a_rr_res(length(rr_mean_t)+1:end) = [];
else
    a_rr_res = [a_rr_res;
                a_rr_res(end)*ones((length(rr_mean_t)-length(a_rr_res)), 1)];
end
low_amp_mask =  a_rr_mean_t <= 0.05 | a_rr_res  <= 0.05;

rr_mean_t = ArtifactCorrection(rr_mean_t, low_amp_mask, t_mean);
a_rr_mean_t(rr_mean_t == 0) = 0;

rr_mean = [start_pad; interp1(t_mean, rr_mean_t, t_res, method); ones(size(tail))];
a_rr_mean = [start_pad; interp1(t_mean, a_rr_mean_t, t_res, method); ones(size(tail))];

%% Fuse RR
d_rr_inst = abs([diff(rr_inst); 0]); 
d_rr_mean = abs([diff(rr_mean); 0]);

rr_inst_jump_ids = [1; find(d_rr_inst > 5); length(rr_inst)];
rr_mean_jump_ids = [1; find(d_rr_mean > 5); length(rr_inst)];
rr_jumps_ids = sort(unique([rr_inst_jump_ids; rr_mean_jump_ids])); 

rr = nan(size(rr_inst));
a_rr = nan(size(rr_inst));
for i=2:length(rr_jumps_ids)
    range = (rr_jumps_ids(i-1)+1):rr_jumps_ids(i);
    mean_diff = abs(diff([median(rr_inst(range)), median(rr_mean(range))]));
%     var_diff = abs(diff([var(rr_inst(range)), var(rr_mean(range))]));
    if mean_diff < 3
        rr(range) = mean([rr_inst(range) rr_mean(range)]');
        a_rr(range) = mean([a_rr_inst(range) a_rr_mean(range)]');
    else
        diffs = rr_jumps_ids(i-1) - rr_mean_jump_ids;
        diffs(diffs < 0) = inf;
        [~, id] = min(diffs);
        mean_sec_start = rr_mean_jump_ids(id);
        diffs = rr_mean_jump_ids - rr_jumps_ids(i);
        diffs(diffs < 0) = inf;
        [~, id] = min(diffs);
        mean_sec_end = rr_mean_jump_ids(id);
        mean_sec_len = diff([mean_sec_start, mean_sec_end])/fs;
        if isempty(mean_sec_len)
            mean_sec_len = 0;
        end

        diffs = rr_jumps_ids(i-1)- rr_inst_jump_ids;
        diffs(diffs < 0) = inf;
        [~, id] = min(diffs);
        inst_sec_start = rr_inst_jump_ids(id);
        diffs = rr_inst_jump_ids - rr_jumps_ids(i);
        diffs(diffs < 0) = inf;
        [~, id] = min(diffs);
        inst_sec_end = rr_inst_jump_ids(id);
        inst_sec_len = diff([inst_sec_start, inst_sec_end])/fs;
        if isempty(inst_sec_len)
            inst_sec_len = 0;
        end

%         if max(inst_sec_len, mean_sec_len) > 10
% %             inst_sec_weigth = (log10((median(rr_inst(inst_sec_start:inst_sec_end))+1)/1)+1);
% %             mean_sec_weigth = (log10((median(rr_mean(mean_sec_start:mean_sec_end))+1)/1)+1);
%             if inst_sec_len > mean_sec_len && median(rr_inst(range)) > 0
%                 rr(range) = rr_inst(range);
%             else
%                 rr(range) = rr_mean(range);
%             end
%         end
        if inst_sec_len > 10 && mean_sec_len > 10
            if median(rr_inst(inst_sec_start:inst_sec_end)) > median(rr_mean(mean_sec_start:mean_sec_end))
                rr(range) = rr_inst(range);
                a_rr(range) = a_rr_inst(range);
            else
                rr(range) = rr_mean(range);
                a_rr(range) = a_rr_mean(range);
            end
        elseif inst_sec_len > 10
            rr(range) = rr_inst(range);
            a_rr(range) = a_rr_inst(range);
        elseif mean_sec_len > 10
            rr(range) = rr_mean(range);
            a_rr(range) = a_rr_mean(range);
        end
    end
end

% Fill missing values, if no to long gap. Otherwise set to zero
rr = fillmissing(rr,'movmedian', 5,'SamplePoints',t_res,'MaxGap',5);
rr(isnan(rr)) = 0;
a_rr(isnan(a_rr)) = 0;
%% compute moving variance
var_rr = movvar(rr, win_length*fs);
var_a_rr = movvar(a_rr, win_length*fs);

%% remove padding and add potential start delay
t = t + t_start - pad_length;

up_sig = [nan(round(t_start*fs*upfac),1); up_sig((pad_length)*fs*upfac:end-pad_length*fs*upfac)];

rr_inst = undoPadAndCrop(rr_inst, fs, pad_length, t_start, sig_len);
a_rr_inst = undoPadAndCrop(a_rr_inst, fs, pad_length, t_start, sig_len);
w_rr = undoPadAndCrop(w_rr, fs, pad_length, t_start, sig_len);

rr_mean = undoPadAndCrop(rr_mean, fs, pad_length, t_start, sig_len);
a_rr_mean = undoPadAndCrop(a_rr_mean, fs, pad_length, t_start, sig_len);

a_rr = undoPadAndCrop(a_rr, fs, pad_length, t_start, sig_len);
rr = undoPadAndCrop(rr, fs, pad_length, t_start, sig_len);
var_rr = undoPadAndCrop(var_rr, fs, pad_length, t_start, sig_len);
var_a_rr = undoPadAndCrop(var_a_rr, fs, pad_length, t_start, sig_len);

%% plot results
if UsePlot
    h = [];
    ax = [];

    % signal and detected temporal characteristics
    h(end+1) = figure();
    findpeaks(up_sig, fs*upfac,'Annotate', 'extents', findpeaks_config{:});
    ax(end+1) = gca;

    % spectral characteristics
    if strcmpi(SpectrogramType, 'fft')
        ShowSpectrogram(spec_sig, spec_fs, spectrogram_config{:}, ...
                        'Type', 'imagesc', ...
                        'Scale', 'logarithmic',...
                        'ColorMap', 'gray', ...
                        'FrequencyUnit','BPM', ...
                        'TimeDelta', t_start - pad_length);
    else
        h_f = figure;
        imagesc(t, f*60, P);
        colormap('pink');

        xlim([t(1), t(end)])
        ylim([min(f(:)), max(f(:))]*60)

        h_c = colorbar();
        h_c.Label.Interpreter = 'latex';
        h_c.Label.String = '$$\frac{power}{frequency} \left[\frac{1}{Hz}\right]$$';

        ax = gca;
        ax.Layer = 'top';
        ax.YDir = 'normal';
        ax.YScale = 'log';

        xlabel('time [s]');
        ylabel('frequency [BPM]');
    end
    ShowSignal(rr_inst, fs, 'PlotParams', 'r', 'NewFigure', false);
    hold on; 
    plot(t_mean, rr_mean_t, 'b');
    h(end+1) = gcf;
    ax(end+1) = gca;

    % RR
    h(end+1) = ShowSignal({rr_inst, rr_mean, rr}, fs);
    ax(end+1) = gca;
    title('Respiratory Rate');

    % resp. amplitude
    h(end+1) = ShowSignal({a_rr_inst, a_rr_mean, a_rr}, fs);
    ax(end+1) = gca;
    title('Amplitude');

    % resp. width
    h(end+1) = ShowSignal(w_rr, fs);
    ax(end+1) = gca;
    title('Width');

    % resp. var
    h(end+1) = ShowSignal({var_rr, var_a_rr}, fs);
    ax(end+1) = gca;
    title('Variances');

    linkaxes(ax, 'x');
    xlim([t(1) t(end)]); 
    
    if AutoClose
        autoArrangeFigures(2, 3, 1)
        res = input('Close figures? [Y]/n:    ','s');
        if ~strcmpi(res,'n')
            for i=1:length(h)
                if ishandle(h(i))
                    close(h(i)); 
                end
            end
        end
    end
end

function [id] = maxID(x)
    [~,id] = max(x);
end

function [sig] = undoPadAndCrop(sig, fs, pad_length, t_start, sig_len)
    sig = [nan(round(t_start*fs),1); sig((pad_length)*fs:end-pad_length*fs)];
    dif_len = sig_len - length(sig);
    if dif_len > 0
        sig = [sig; nan(dif_len, 1)];
    else
        sig = sig(1:sig_len);
    end
end

function [rr] = ArtifactCorrection(rr, low_amp_mask, t)    
    % Filter values, which jump to much
    d_rr = [0; diff(rr)]./[0; diff(t)];  
    jump_ids = [1; find(abs(d_rr) > 5); length(rr)];
    mask = false(length(rr), 1);
    for i=2:length(jump_ids)
        if (t(jump_ids(i)) < t(jump_ids(i-1)) + jump_eps)
            mask(jump_ids(i-1):jump_ids(i)) = 1;
        end
    end
    rr(mask) = nan;
    
    % Filter values with to little amplitude, if med is different or the variance is to high
    low_amp_mask(isnan(rr)) = 0;
    ids = [1; find(abs(diff(low_amp_mask))>0)+1; length(rr)];
    amp_low = low_amp_mask(ids);
    for i=1:length(ids)-1
        range = ids(i):(ids(i+1)-1);
        if amp_low(i)
            if isempty(range)
                range = ids(i);
            end
            sec_vals = rr(range);
            sec_vals(isnan(sec_vals)) = [];
            if isempty(sec_vals)
                continue;
            end
            sec_med = nanmedian(sec_vals);
            sec_var = diff([min(sec_vals), max(sec_vals)])/sec_med;
            confidence_len = 5*60/sec_med;
            if  sec_var > 0.45 || sec_med < 5 || diff([t(range(1)), t(range(end))]) < confidence_len
                rr(range) = nan;
            end
        end
    end
    
    % filter values, which are only short segments between 0 or nan
    zero_nan_mask = (isnan(rr) | rr == 0)';
    ids = [find(abs(diff(zero_nan_mask))>0) length(rr)];
    zero_or_nan = ~zero_nan_mask(ids);
    for i=1:length(ids)-1
        range = [ids(i), (ids(i+1))];
        if ~zero_or_nan(i)
            sec_len = diff(t(range));
            if  sec_len < jump_eps
                rr(range(1):range(2)) = nan;
            end
        end
    end
    
    % Fill missing values, if no to long gap. Otherwise set to zero
    rr = fillmissing(rr,'movmedian', 5,'SamplePoints',t,'MaxGap',10);
    rr(isnan(rr)) = 0;
end
end


