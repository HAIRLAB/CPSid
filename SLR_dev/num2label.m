function label = num2label(label_binary, label_names)
% REPLACE A BINARY LABEL WITH LABEL_NAMES
%
% -- Usage
% label = num2label(label_binary, label_names);
%
% -- Input
% label_binary : label vector consisting of [1 2]
%
% 2013/11/1 O.Yamashita
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.


if iscell(label_names)
    label = cell(length(label_binary),1);
    ix = find(label_binary == 1);
    for ii = 1 : length(ix)
        label{ix(ii)} = label_names{1};
    end
    ix = find(label_binary == 2);
    for ii = 1 : length(ix)
        label{ix(ii)} = label_names{2};
    end

elseif isnumeric(label_names)
   
    label = zeros(length(label_binary),1);
    label(label_binary==1) = label_names(1);
    label(label_binary==2) = label_names(2);
    
else
    error(' Class of ''label_names'' is either string or numeric !');
    
end