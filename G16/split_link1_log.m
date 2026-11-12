function files = split_link1_log(logfile, opts)
%SPLIT_LINK1_LOG  Split a Gaussian log of chained (--Link1--) jobs into one log per job.
%
%   split_link1_log('NMA_WB97XD_DEF2TZVP_TD_ALL.LOG')
%   split_link1_log('NMA_wB97XD_def2TZVP_scan_all.log')        % field scans
%   files = split_link1_log(..., 'naming', 'chk')
%
%   Two ways of naming the output logs ('naming' option):
%     'auto'  (default) field-scan title if found in the first job
%             ("... scan along X: N = ...", as written by make_field_scan.m),
%             otherwise the %chk name of each job.
%     'field' <prefix>_Fp020.log, <prefix>_Fm040.log, ... from the title
%             (old behaviour; error if the title is missing).
%     'chk'   <name of the job's .chk file>.log, e.g.
%             NMA_wB97XD_def2TZVP_TD_trans.log, NMA_wB97XD_IP_EA ... .log
%
%   The internal Link1 step that Gaussian creates for the Freq part of an
%   Opt+Freq job carries the same title / %chk (or none) and is kept
%   together with its optimization.
%
%   Name-value options:
%     naming  : 'auto' | 'field' | 'chk'          (default 'auto')
%     pattern : regexp with one token capturing N (default 'scan along X: N = (-?\d+)')
%     prefix  : field naming only; default = log name without trailing '_scan_all'
%     outdir  : default = folder of the input log
%
%   Output: cell array with the full names of the logs written.

arguments
    logfile (1,:) char
    opts.naming (1,:) char {mustBeMember(opts.naming, {'auto','field','chk'})} = 'auto'
    opts.pattern (1,:) char = 'scan along X: N = (-?\d+)'
    opts.prefix (1,:) char = ''
    opts.outdir (1,:) char = ''
end

[p, b] = fileparts(logfile);
if isempty(opts.outdir)
    if isempty(p), opts.outdir = pwd; else, opts.outdir = p; end
end
if isempty(opts.prefix)
    opts.prefix = regexprep(b, '_scan_all$', '', 'ignorecase');
end
if ~exist(opts.outdir, 'dir'), mkdir(opts.outdir); end

lines = splitlines(fileread(logfile));
L = numel(lines);

% ---------- segments between "Entering Link 1" lines ----------
starts = find(contains(lines, 'Entering Link 1'));
if isempty(starts)
    error('No "Entering Link 1" line found in %s.', logfile);
end
starts(1) = 1;                      % keep the file header with the first job
ends = [starts(2:end) - 1; L];
nseg = numel(starts);

% ---------- key of each segment: field value N or %chk name ----------
segN   = nan(nseg, 1);
segChk = strings(nseg, 1);
for s = 1:nseg
    blk = lines(starts(s):ends(s));
    tok = regexp(blk, opts.pattern, 'tokens', 'once');
    hit = find(~cellfun(@isempty, tok), 1);
    if ~isempty(hit), segN(s) = str2double(tok{hit}{1}); end
    tok = regexp(blk, '^\s*%chk\s*=\s*(\S+)', 'tokens', 'once', 'ignorecase');
    hit = find(~cellfun(@isempty, tok), 1);
    if ~isempty(hit)
        % Gaussian wraps long Link0 lines at 80 columns: join the continuation
        % lines until the name is complete (ends in .chk)
        cname = strtrim(tok{hit}{1});
        j = hit + 1;
        while isempty(regexpi(cname, '\.chk$', 'once')) && j <= numel(blk) && j <= hit + 4
            nxt = strtrim(blk{j});
            if isempty(nxt) || startsWith(nxt, {'-', '%', '#'}), break; end
            cname = [cname nxt]; %#ok<AGROW>
            j = j + 1;
        end
        cname = strrep(strrep(cname, '\', '/'), '/', filesep);
        [~, cb] = fileparts(cname);
        segChk(s) = cb;
    end
end

naming = opts.naming;
if strcmp(naming, 'auto')
    if ~isnan(segN(1)), naming = 'field'; else, naming = 'chk'; end
end
switch naming
    case 'field'
        if isnan(segN(1))
            error(['Title pattern not found in the first job. Check ''pattern'', ' ...
                   'or use ''naming'',''chk''.']);
        end
        for s = 2:nseg, if isnan(segN(s)), segN(s) = segN(s-1); end, end
        key = string(segN);
    case 'chk'
        if segChk(1) == ""
            error(['No chk line (%chk=...) found in the first job: use ''naming'',''field'' ' ...
                   'or add a %chk line to the input.']);
        end
        for s = 2:nseg, if segChk(s) == "", segChk(s) = segChk(s-1); end, end
        key = segChk;
end

% ---------- group consecutive segments with the same key ----------
gStart = [1; find(key(2:end) ~= key(1:end-1)) + 1];
gEnd   = [gStart(2:end) - 1; nseg];

fprintf('%-46s %5s %8s %s\n', 'File', 'Job', 'Normal', 'Status');
written = {};
files = cell(numel(gStart), 1);
for g = 1:numel(gStart)
    a = starts(gStart(g));  z = ends(gEnd(g));
    if strcmp(naming, 'field')
        N = segN(gStart(g));
        if N >= 0, name = sprintf('%s_Fp%03d', opts.prefix, N);
        else,      name = sprintf('%s_Fm%03d', opts.prefix, -N); end
    else
        name = char(key(gStart(g)));
    end
    base = name;  k = 1;            % avoid overwriting if a key appears twice
    while any(strcmp(written, name)), k = k + 1; name = sprintf('%s_%d', base, k); end
    written{end+1} = name;          %#ok<AGROW>

    if strcmpi(fullfile(opts.outdir, [name '.log']), fullfile(p, [b '.log']))
        name = [name '_job'];       % never overwrite the input log
    end
    files{g} = fullfile(opts.outdir, [name '.log']);
    fid = fopen(files{g}, 'wt');
    fprintf(fid, '%s\n', lines{a:z});
    fclose(fid);

    blk   = lines(a:z);
    nNorm = sum(contains(blk, 'Normal termination'));
    nLink = gEnd(g) - gStart(g) + 1;        % Gaussian link steps in this job
    if any(contains(blk, 'Error termination'))
        status = 'ERROR termination';
    elseif nNorm >= nLink
        status = 'OK';
    elseif nNorm > 0
        status = sprintf('%d of %d steps terminated', nNorm, nLink);
    else
        status = 'incomplete / still running';
    end
    if any(contains(blk, 'Stationary point found')), status = [status ', stationary point']; end
    fprintf('%-46s %5d %8d %s\n', [name '.log'], g, nNorm, status);
end
fprintf('%d job logs written to %s\n', numel(gStart), opts.outdir);
end
