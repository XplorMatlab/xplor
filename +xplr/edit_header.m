function header = edit_header(dat)
% function header = edit_header(dat)
%---
% Open graphic interface (headerEdit object) to let user define header
% information for ND-array dat.

% graphic interface for setting headers
E = xplr.HeaderEdit(dat);

% wait for 'ok' button press (i.e. when the interface figure closes)
waitfor(E.hf)

% get headers
header = E.header;
