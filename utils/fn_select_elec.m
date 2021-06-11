function elec = fn_select_elec(cfg, elec_og)

% SELECT_ELEC selects a subset of electrodes based on their labels
%
% Configuration options:
%   cfg.channel     = Nx1 cell array of channels to select


% get the default options
channel = ft_getopt(cfg, 'channel');
chansel = ft_channelselection(channel, elec_og.label);

if isstr(channel) % if channel is given as a string for one electrode, change it to a cell
  channel = {channel};
end

% check to make sure that the target channels exist in the elec structure
for n = 1:length(channel)
  if max(unique(isstrprop(channel{n}, 'digit'))); % only check if the target channels refer to specifically numbered electrodes
    if ~any(match_str(chansel, channel{n}))
      warning('Some of the target channels do not exist in the original electrode structure, so select_elec is only taking those that do exist')
    end
  end
end

% OLD WORKAROUND for making function compatible with specific electrode
% inputs with numbers or with wildcards for numbers to select all
% electrodes in a grid or depth tract
%
% for n = 1:numel(chansel)
%   if any(cell2mat(strfind(chansel, '1'))) || any(cell2mat(strfind(chansel, '2'))) || any(cell2mat(strfind(chansel, '3'))) || ...
%       any(cell2mat(strfind(chansel, '4'))) || any(cell2mat(strfind(chansel, '5'))) || any(cell2mat(strfind(chansel, '6'))) || ...
%       any(cell2mat(strfind(chansel, '7'))) || any(cell2mat(strfind(chansel, '8'))) || any(cell2mat(strfind(chansel, '9')))
%     % do nothing
%   else
%     tmp = chansel{n};
%     chansel{n} = [];
%     for m = 1:100;
%       chansel{end+1} = [tmp num2str(m)];
%     end
%   end
% end

elec = [];

% Pre-allocate the subfields in  elec accordingly
sub_vals = fieldnames(elec_og);
for n = 1:numel(sub_vals)
  if eval(['ischar(elec_og.' sub_vals{n} ')'])
    eval(['elec.' sub_vals{n} '= elec_og.' sub_vals{n} ';']);
  elseif eval(['iscell(elec_og.' sub_vals{n} ')'])
    eval(['elec.' sub_vals{n} '= {};']);
  elseif eval(['ismatrix(elec_og.' sub_vals{n} ')'])
    eval(['elec.' sub_vals{n} '= [];']);
  elseif eval(['isstruct(elec_og.' sub_vals{n} ')'])
    eval(['elec.' sub_vals{n} '= elec_og.' sub_vals{n} ';']);
  end
end

for e = 1:numel(elec_og.label) % for each electrode in original structure
  for s = 1:numel(chansel) % for each target electrode
    if match_str(char(elec_og.label{e}), chansel{s}) % if the target electrode is in the original structure
      for n = 1:numel(sub_vals) % for each of the subfields in the elec structure
        if eval(['size(elec_og.' sub_vals{n} ',1) == numel(elec_og.label)']) && numel(elec_og.label) ~=1 % if this subfield describes something that looks like it describes electrodes on an elec by elec basis
          if eval(['~isempty(elec_og.' sub_vals{n} ') && iscell(elec_og.' sub_vals{n} ')'])
            eval(['elec.' sub_vals{n} '(end+1,1) = elec_og.' sub_vals{n} '(e);'])
          elseif eval(['~isempty(elec_og.' sub_vals{n} ') && (ismatrix(elec_og.' sub_vals{n} ') && ~isstruct(elec_og.' sub_vals{n} ') && ~ischar(elec_og.' sub_vals{n} '))'])
            eval(['elec.' sub_vals{n} '(end+1, :) = elec_og.' sub_vals{n} '(e, :);'])
          end
        end
      end
    end
  end
end

% copy over fields that do not appear to be related to each electrode
for n = 1:numel(sub_vals)
  if eval(['size(elec_og.' sub_vals{n} ',1) ~= numel(elec_og.label)'])
    eval(['elec.' sub_vals{n} '= elec_og.' sub_vals{n} ';']);
  end
end

% for depths and strips
for n = 1:numel(sub_vals) % fieldnames loop
  if strcmp(sub_vals{n}, 'surface') || strcmp(sub_vals{n}, 'depth');
    % loop through each subset and remove the labels not in cfg.channel
    for s = 1:numel(eval(['elec.' sub_vals{n}])) % subsetted grid/strip/depth loop
      tmp = eval(['elec.' sub_vals{n}]);
      subset = tmp{s};
      
      for e = numel(subset):-1:1 % electrode loop for electrodes within this subset
        if isempty(match_str(elec.label, subset{e}))
          subset(e) = [];
        end
      end
      
      tmp{s} = subset;
      
      eval(['elec.' sub_vals{n} '= tmp;'])
    end
  end
end

% remove the tra field if it exists
if isfield(elec_og, 'tra')
  elec = rmfield(elec, 'tra');
end

