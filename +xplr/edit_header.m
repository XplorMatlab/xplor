function header = edit_header(dat)
% function header = edit_header(dat)
%---
% Open graphic brick.interface (H   eaderEdit object) to let user define header
% information for ND-array dat.

% graphic brick.interface for setting headers
E = xplr.HeaderEdit(dat);

% wait for 'ok' button press (i.e. when the brick.interface figure closes)
waitfor(E.hf)

% get headers
header = E.header;
