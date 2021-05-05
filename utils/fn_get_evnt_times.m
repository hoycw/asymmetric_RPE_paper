function [event_times] = fn_get_evnt_times(an_event_type,plt_event_lab)
%% Find the time (sec) of an event relative to the time-locked events for plotting
% INPUTS:
%   an_event_type [str] - label of the time-locked event
%       {'S','R','F'} for stimulus, response, feedback (FB)
%   plt_event_lab [cell array] - set of labels of events to plot
%       {'S','R','F','Fon','Foff'} for stimulus, response, feedback, FB onset, FB offset
% OUTPUTS:
%   event_times [array] - times of the events (in sec) relative to analysis time-locked event

event_times = zeros(size(plt_event_lab));
if strcmp(an_event_type,'S')
    for evnt_ix = 1:numel(plt_event_lab)
        switch plt_event_lab{evnt_ix}
            case 'S'
                event_times(evnt_ix) = 0;
            case 'R'
                error('this looks like a bug, would need to load and add bhv.rt');
                event_times(evnt_ix) = bhv.prdm.target;
            case {'F','Fon'}
                event_times(evnt_ix) = bhv.prdm.target+bhv.prdm.fb_delay;
            case 'Foff'
                event_times(evnt_ix) = bhv.prdm.target+bhv.prdm.fb_delay+bhv.prdm.fb;
            otherwise
                error(['Unknown event type in plt: ' plt_event_lab{evnt_ix}]);
        end
    end
elseif any(strcmp(an_event_type,{'F','R'}))
    event_times(1) = 0;
else
    error('Unknown an_event_type');
end

end