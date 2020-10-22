function host = hostname()
%HOSTNAME Return an identifiant specific to the computer in use
%---
% function hostname = hostname()
%---
% returns an identifiant specific to the computer in use

% Thomas Deneux
% Copyright 2015-2017

comp = computer;
switch comp
    case {'PCWIN' 'PCWIN64'}
        comp = 'PCWIN';
        host = getenv('COMPUTERNAME');
    otherwise
        host = getenv('HOSTNAME');
        if isempty(host)
            [dum host] = system('echo $HOSTNAME');
            host = strrep(host,char(10),''); % remove endlines
        end %#ok<*ASGLU>
end
host = [comp '-' host];

