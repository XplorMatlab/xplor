function header = editHeader(dat)
% function header = editHeader(dat)
%---
% Open graphic interface (headerEdit object) to let user define header
% information for ND-array dat.

% graphic interface for setting headers
E = xplr.headerEdit(dat);

% wait for 'ok' button press (i.e. when the interface figure closes)
waitfor(E.hf)

% get headers
header = E.header;

% register them to the bank
if ~isempty(header)
    xplr.bank.registerheaders(header)
end


