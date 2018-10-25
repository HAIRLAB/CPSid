function [nx, scale, base] = normalize_feature(x, scale_mode, mean_mode, scale, base, verbose)
% NORMALIZE FEATURE VECTORS BEFORE CLASSIFICATION.
% First mean adjustment, then scaling. 
% 
% -- Usage
% [nx, scale, base] = normalize_feature(x, scale_mode, mean_mode)
% [nx, scale, base] = normalize_feature(x, scale_mode, mean_mode,
% scale, base)
%
% -- Input
% x : input features ( trial * feature )
% scale_mode : {all', 'each', 'stdall', 'stdeach', 'none', 'medeach' }
% mean_mode : {all', 'each', 'none', 'medeach' }
% (scale, base) : optional
%
% -- Example
% > [nx] = normalize_feature(x, 'stdeach', 'each') 
% Each column in 'x' is normalized to have mean 0 and std 1.
%
% 2011/06/25 Okito Yamashita
% * add the case 'scale' or 'base' is empty 
% 2011/04/08 Okito Yamashita
% * add options to robust base and scale estimates to outliers
% * add "verbose" not to display comments
% 2007/04/05 Okito Yamashita
% * 'stdall' and 'stdeach' are added to 'scale_mode' 
% 2006/01/31 Okito Yamashita
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

if nargin < 6,
    verbose = 0;
end
    

%% mean adjustment
switch mean_mode     
    case 'all',
        if nargin < 5 | isempty(base)
        base = mean(x(:));
        end
        nx = x - base;
        if verbose, fprintf(' Mean adjustment by all features mean !'); end
    case 'each',
        if nargin < 5 | isempty(base)
        base =  mean(x,1);  % adjust baseline
        end
        nx = x - repmat(base, [size(x,1),1]);
        if verbose, fprintf(' Mean adjustment by each feature mean!'); end
    case 'medeach',
        if nargin < 5 | isempty(base)
        base =  median(x,1);  % adjust baseline
        end
        nx = x - repmat(base, [size(x,1),1]);
        if verbose, fprintf(' Mean adjustment by each feature median!'); end
    otherwise,
        if nargin < 5 | isempty(base)
        base = 0;
        end
        nx = x;
   %     fprintf(' No mean adjustment !');
end

%% scale adjustment
switch scale_mode     
    case 'all',   % scaling by abs(max)
        if nargin < 4 | isempty(scale)
        scale = max(abs(nx(:))); 
        end
        nx = nx / scale ;
        if verbose, fprintf(' Scaling by all feature maximum !\n'); end
    case 'each',    % scaling by abs(max)
        if nargin < 4 | isempty(scale)
        scale = max(abs(nx),[],1); 
        end
        nx = nx ./repmat(scale, [size(nx,1),1]) ;
        if verbose, fprintf(' Scaling by each featue maximum !\n'); end
    case 'stdall',   % scaling by std.
        if nargin < 4 | isempty(scale)
        scale = std(nx(:),[],1);  % scalar 
        end
        nx = nx / scale ;
        if verbose, fprintf(' Scaling feature by all feature standard deviation !\n'); end
    case 'stdeach', % scaling by std.
        if nargin < 4 | isempty(scale)
        scale = std(nx,[],1); 
        end
        nx = nx ./repmat(scale, [size(nx,1),1]) ;
        if verbose, fprintf(' Scaling feature by featue by its standard deviation !\n'); end
    case 'medeach'
        if nargin < 4 | isempty(scale)
            scale = 1.4826 * median(abs(nx-repmat(median(nx,1),[size(nx,1),1])),1);  % median absolute deviation (robust deviation estimation)
        end
        nx = nx ./repmat(scale, [size(nx,1),1]) ;
        if verbose, fprintf(' Scaling feature by featue by its median absolute deviation !\n'); end
        
    otherwise,
        if nargin < 4 | isempty(scale)
        scale = 1;
        end
        nx = nx;
    %    fprintf(' No scaling !\n');
end



