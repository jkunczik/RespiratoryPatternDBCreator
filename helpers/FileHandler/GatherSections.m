function [sections] = GatherSections(sigs)
%GATHERSECTIONS 
    
    labels = enumeration('Pattern'); 

    for i=1:length(labels)-1
        secs = [];
        for j=1:length(sigs)
            new_secs = GetSection(sigs(j), labels(i));
            if (labels(i) == Pattern.BIOT && length(new_secs) < 2)
                % because both BIOT secstions are close to each other (only
                % 5 second invaldi between them), they might accidentally 
                % be interpreted as one
                new_secs = [new_secs, new_secs];
                middle = round(length(new_secs(1).sig)/2);
                new_secs(1).sig = new_secs(1).sig(1:middle,:);
                new_secs(2).sig = new_secs(2).sig((middle+1):end,:);
                new_secs(2).start_id = new_secs(2).start_id + middle + 1;
            end
            secs = [secs, new_secs];
        end
        eval(sprintf('sections.%s = secs;', labels(i)));
    end

end

