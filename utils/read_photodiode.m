function [event_onsets, event_offsets, data_shades] = read_photodiode(data_photo_orig, min_event_length, n_shades, plot_it)
% This function takes an array of raw photodiode values and outputs event onsets and offsets
% 
% NOTE: The lowest values are considered non-events. This function assumes positive values for ON events.
% 
% data_photo_orig  - [1xN] array of photodiode values. This code assumes that photodiode events are
%                     positive values (white). If you use black for photodiode events, then multiply
%                     your photodiode data by -1 before using this function.
% min_event_length - scalar. This is in points, not seconds. Any changes in the photodiode data with a
%                     duration less than this will be ignored as noise.
% n_shades         - scalar. The number of possible shades, or event types. If you are only using
%                     black and white, set this to zero.
% 
% event_onsets   - A cell array of arrays. The index of the cell array is the shade number. The index
%                   of the interior arrays corresponds to the event number. The contents of each array
%                   is time points of the events.
% event_offsets  - Same as for event_onsets, but for event offsets.
% data_shades    - [1xN] array of values corresponding to the raw photodiode data. The value of each
%                   time point is the event number.
% plot_it        - 0/1 flag for plotting



n_hist_bins = 100; % Set this to a higher number for more shades. More bins can mean more errors though.

% Transpose if necessary
if size(data_photo_orig,1) > size(data_photo_orig,2)
  data_photo_orig = data_photo_orig';
end

% Dilate and then erode the mic amplitudes to reduce noise
%  Dilate is kind of like doing the moving-maximum. Eroding then makes it "zero-phase".
%   This is also know as the morphological open in image processing.
%   These functions require the image processing toolbox
filt_window = ones(1,floor(min_event_length*0.500));
data_photo = imdilate(data_photo_orig,filt_window);
data_photo = imerode(data_photo,filt_window);

% Bin the data amplitudes to determine cutoff points between different shades
[hist_counts, hist_edges] = hist(data_photo,n_hist_bins);
% Add zeros to edges to allow for peaks at edges to be found
hist_counts = [0 hist_counts 0];
hist_edges = [hist_edges(1)-(hist_edges(2)-hist_edges(1)) hist_edges hist_edges(end)+(hist_edges(2)-hist_edges(1))];

% Find the top n_shades peaks in the histogram. These will correspond to the amplitude levels for the different shades
[~, peak_idx] = findpeaks(hist_counts,'MinPeakDistance',min(n_hist_bins/20),'SortStr','descend');
top_n = sort(hist_edges(peak_idx(1:n_shades)),'ascend');

% Set the minimum and maximum amplitude for each shade
%   shade number one is the lowest value, so no cutoff
cutoff_vals = zeros(1,n_shades);
for shade_n = 2:n_shades
  cutoff_vals(shade_n) = top_n(shade_n-1) + (top_n(shade_n)-top_n(shade_n-1))/4; % Dividing by 4 to be closer to lower edge
end

% Set event data to detected shade numbers
data_shades = zeros(size(data_photo));
for shade_n = 1:n_shades
  % Increment value in data_shades if data_photo at that time point is greater than threshold
  data_shades(data_photo >= cutoff_vals(shade_n)) = data_shades(data_photo >= cutoff_vals(shade_n)) + 1;
end

% Calculate first derivative to get edges which correspond to onsets and offsets
data_events = [diff(data_shades) 0]; % Add a point because diff removes one

% If an edge is not sharp enough, this can lead to an intermediate shade being assigned near the real shade.
%   Check for this and correct
for sample_n = 1:(length(data_events)-min_event_length)
  if ~isempty(find((data_events(sample_n:(sample_n+min_event_length))),1))
    % There is at least one event in this time range
    if length(find((data_events(sample_n:(sample_n+min_event_length))))) > 1
      data_tmp = data_events(sample_n:(sample_n+min_event_length));
      if sum(data_tmp) > 0
        % Positive slope. Push maximum value towards the earlier time period
        first_nz = find(data_tmp,1);
        data_sum = sum(data_tmp);
        data_tmp = zeros(size(data_tmp));
        data_tmp(first_nz) = data_sum;
      end
      if sum(data_tmp) < 0
        % Negative slope. Push maximum value towards the later time period
        last_nz = find(data_tmp,1,'last');
        data_sum = sum(data_tmp);
        data_tmp = zeros(size(data_tmp));
        data_tmp(last_nz) = data_sum;
      end
      % Mixed positive and negative should not be possible from morphological open
      data_events(sample_n:(sample_n+min_event_length)) = data_tmp;
    end
  end
end

% Recreate the data_shades array
data_shades = data_shades(1) + cumsum(data_events); % Adding first value of original because diff removes constant

% Put onsets and offsets into cell arrays
event_onsets{1} = []; % Zero is considered no event
for shade_n = 2:n_shades
  event_onsets{shade_n} = find(data_events==shade_n-1)';
end
event_offsets{1} = []; % Zero is considered no event
for shade_n = 2:n_shades
  event_offsets{shade_n} = find(-data_events==shade_n-1)';
end

% Plot data, histogram of photo amplitudes, detected peaks, cutoff values, etc
if plot_it
    figure;
    subplot(3,1,1);
    plot(data_photo_orig,'k'); hold on;
    for shade_n = 2:n_shades
        for event_n = 1:length(event_onsets{shade_n})
            plot([event_onsets{shade_n}(event_n) event_onsets{shade_n}(event_n)],[min(data_photo_orig) max(data_photo_orig)],'g');
        end
        for event_n = 1:length(event_offsets{shade_n})
            plot([event_offsets{shade_n}(event_n) event_offsets{shade_n}(event_n)],[min(data_photo_orig) max(data_photo_orig)],'r');
        end
    end
    title('Raw photodiode data with event onsets');
    subplot(3,1,2);
    plot(data_shades, 'k'); ylim([0.5 n_shades+0.5]);
    title('Detected shades');
    subplot(3,1,3);
    bar(hist_edges,hist_counts); hold on;
    plot(hist_edges(peak_idx(1:n_shades)), hist_counts(peak_idx(1:n_shades)), 'r*');
    for shade_n = 2:n_shades
        plot([cutoff_vals(shade_n) cutoff_vals(shade_n)],[0 hist_counts(peak_idx(1))],'r');
    end
    title('Histogram of amplitude values, peaks, and cutoffs');
end


