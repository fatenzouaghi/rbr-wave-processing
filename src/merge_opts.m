function out = merge_opts(defaults, user)
% MERGE_OPTS  Struct overlay: user values override defaults.
out = defaults;
if isempty(user), return; end
fn = fieldnames(user);
for i = 1:numel(fn)
    out.(fn{i}) = user.(fn{i});
end
end
