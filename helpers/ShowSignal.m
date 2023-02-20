function h = ShowSignal( signal, fs, varargin )
%SHOWSIGNAL Plots a continously sampled signal
%   if signal is a cell array of signals, the function plots all signals
%   recursively

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('NewFigure', true);
p.addParameter('PlotParams','');
p.addParameter('SimilarSignals',true);
p.addParameter('Label', 'amplitude [a. u.]');
p.addParameter('FullScreen', false);
p.addParameter('SubPlot',[]);
p.parse(varargin{:});

NewFigure = p.Results.NewFigure;
PlotParams = p.Results.PlotParams;
SimilarSignals = p.Results.SimilarSignals;
Label = p.Results.Label;
SubPlot = p.Results.SubPlot;

if iscell(signal)
    if NewFigure
        figs = figure; 
    end
    for i=1:length(signal)
        if ~SimilarSignals
            subplot_opts = [length(signal) 1 i];
        else
            subplot_opts = [];
        end
        ShowSignal( signal{i}, fs, 'NewFigure', false,...
                                   'PlotParams', PlotParams,...
                                   'SimilarSignals',true,...
                                   'Label', Label,...
                                   'SubPlot', subplot_opts);
    end
    if nargout >= 1
        h = figs; 
    end
    return;
end

if NewFigure
    figs = figure; 
end

sigSize = size(signal);
if isempty(sigSize(sigSize == 1))
    %signal array with multiple column signals

    for i = 1:sigSize(2)
        if ~SimilarSignals
            subplot(sigSize(2),1,i);
        end
        ShowSignal( signal(:,i), fs, 'NewFigure', false,...
                                     'PlotParams', PlotParams );
    end
    if nargout >= 1
        h = figs; 
    end
else
    hold on
    if SubPlot
       subplot(SubPlot(1), SubPlot(2), SubPlot(3)); 
    end
    if iscell(PlotParams)
        plot(0:1/fs:(length(signal)-1)*1/fs,signal,PlotParams{:});
    else
        plot(0:1/fs:(length(signal)-1)*1/fs,signal,PlotParams); 
    end
    xlabel('Time [s]', 'Interpreter', 'Latex');
    ylabel(Label, 'Interpreter', 'Latex');
    hold off
end

if p.Results.FullScreen
    figs.Units = 'normalized';
    figs.OuterPosition = [0 0 1 1];
end

if nargout >= 1
   h = figs; 
end

end

