function [nn,ctrs] = hist4(varargin)


[~,args,nargs] = axescheck(varargin{:});

clear cax

if nargs < 1
    error(message('stats:hist4:TooFewInputs'))
end
x = args{1};
 
% See if nbins/ctrs was given as the second argument, or only name/value
% pairs.
if nargs > 1 && ~ischar(args{2})
    binSpec = args{2};
    args = args(3:end); % strip off x and nbins/ctrs
else
    binSpec = [];
    args = args(2:end); % strip off x
end

% Process input parameter name/value pairs, assume unrecognized ones are
% graphics properties.
pnames = {'nbins','ctrs','edges'};
dflts =  { [],     [],       []};
[errid,errmsg,nbins,ctrs,edges] = internal.stats.getargs(pnames, dflts, args{:});
if ~isempty(errmsg)
    error(['stats:hist4:' errid], errmsg);
end

% Make sure they haven't mixed 'nbins'/'ctrs'/'edges' name/value pairs with
% the CTRS or NBINS unnamed second arg syntax, or used more than one of
% those parameter name/value pairs.
if (isempty(nbins)+isempty(ctrs)+isempty(edges)+isempty(binSpec)) < 3
    error(message('stats:hist4:AmbiguousBinSpec'));
elseif ~isempty(binSpec)
    if iscell(binSpec)  % hist4(x,ctrs)
        ctrs = binSpec;
    else                % hist4(x,nbins)
        nbins = binSpec;
    end
end

if ~isempty(nbins)
    % Use the specified number of bars in each direction, centers and edges
    % to be determined.
    histBehavior = true;
    if ~(isnumeric(nbins) && numel(nbins)==3)
        error(message('stats:hist4:BadNbins'));
    end
    autobins = true;
    
elseif ~isempty(ctrs)
    % Use the specified bin centers.
    histBehavior = true;
    if ~(iscell(ctrs) && numel(ctrs)==2 && isnumeric(ctrs{1}) && isnumeric(ctrs{2}))
        error(message('stats:hist4:BadCtrs'));
    end
    ctrs = {ctrs{1}(:)' ctrs{2}(:)' ctrs{3}(:)'};
    autobins = false;
    nbins = [length(ctrs{1}) length(ctrs{2}) length(ctrs{3})];
    
elseif ~isempty(edges)
    % Use the specified bin edges.
    histBehavior = false;
    if ~(iscell(edges) && numel(edges)==2 && isnumeric(edges{1}) && isnumeric(edges{2}))
        error(message('stats:hist4:BadEdges'));
    end
    edges = {edges{1}(:)' edges{2}(:)' edges{3}(:)'};
    autobins = false;
    % Just as with histc, there will be #edges bins
    nbins = [length(edges{1}) length(edges{2}) length(edges{3})];
    
else
    % Assume a 10x10x10 grid of bars, centers and edges to be determined.
    histBehavior = true;
    autobins = true;
    nbins = [10 10 10];
end

[nrows,ncols] = size(x);
if ncols ~= 3
    error(message('stats:hist4:WrongNumCols'));
end

% Special case for empty data (follows what HIST does).
if isempty(x)
    if autobins
       ctrs = {1:nbins(1) 1:nbins(2) 1:nbins(3)};
    end
    %n = zeros(nbins); % Nothing to count, return nbins(1) by nbins(2) zeros
 
else
    % Bin each observation in the x-direction, and in the y-direction.
    bin = zeros(nrows,3);
    for i = 1:3
        minx = min(x(:,i));
        maxx = max(x(:,i));
        
        % If only the number of bins was given, compute edges and centers
        % for equal-sized bins spanning the data.
        if autobins
            if isinf(minx) || isinf(maxx)
                error(message('stats:hist4:InfData'));
            elseif minx == maxx
                minx = minx - floor(nbins(i)/2) - 0.5;
                maxx = maxx + ceil(nbins(i)/2) - 0.5;
            end
            binwidth{i} = (maxx - minx) / nbins(i);
            edges{i} = minx + binwidth{i}*(0:nbins(i));
            ctrs{i} = edges{i}(1:nbins(i)) + binwidth{i}/2;
            % Make histc mimic hist behavior:  everything < ctrs(1) gets
            % counted in first bin, everything > ctrs(end) gets counted in
            % last bin.  ctrs, edges, and binwidth do not reflect that, but
            % histcEdges does.
            histcEdges = [-Inf edges{i}(2:end-1) Inf];
            
        % If the bin centers were given, compute their edges and widths.
        elseif histBehavior
            c = ctrs{i};
            dc = diff(c);
            edges{i} = [c(1) c] + [-dc(1) dc dc(end)]/2;
            binwidth{i} = diff(edges{i});
            % Make histc mimic hist behavior:  everything < ctrs(1) gets
            % counted in first bin, everything > ctrs(end) gets counted in
            % last bin.  ctrs, edges, and binwidth do not reflect that, but
            % histcEdges does.
            histcEdges = [-Inf edges{i}(2:end-1) Inf];
            
        % If the bin edges were given, compute their widths and centers (if
        % asked for).
        else % if ~histBehavior
            e = edges{i};
            de = diff(e);
            histcEdges = e;
            % Make the display mimic bar's histc behavior: an implied bin
            % above edges(end), the same width as the last explicit one.
            % ctrs, edges, and binwidth need that explicitly, histcEdges
            % doesn't.
            edges{i} = [e e(end)+de(end)];
            binwidth{i} = [de de(end)];
            if nargout > 1
                c = zeros(size(de));
                c(1) = e(1) + de(1)/2;
                for j = 2:length(c)
                    c(j) = 2*e(j) - c(j-1);
                end
                % When edges are specified, it may not be possible to return
                % centers for which the edges are midpoints.  Warn if that's
                % the case.
                if any(c <= e(1:end-1)) || ...
                   abs(c(end) - (e(end)-de(end)/2)) > 1000*eps(de(end));
                    warning(message('stats:hist4:InconsistentEdges'));
                    c = e(1:end-1) + de/2;
                end
                ctrs{i} = [c e(end)+de(end)/2];
            end
        end
        
        % Get the 1D bin numbers for this column of x.  Make sure +Inf
        % goes into the nth bin, not the (n+1)th.
        [~,bin(:,i)] = histc(x(:,i),histcEdges,1);
        bin(:,i) = min(bin(:,i),nbins(i));
    end
end
    
% Combine the two vectors of 1D bin counts into a grid of 2D bin
% counts.
 n = accumarray(bin(all(bin>0,2),:),1,nbins);


if 0 < nargout
    nn = n;
    return
end