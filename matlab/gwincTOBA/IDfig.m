function tt = IDfig(varargin)
% IDfig puts the name of the calling function and the date on the right side of a figure
% there must be a current axis, and it must be called from an m-file
% replaces IDplot, which has a name conflict with something in the sys-ID
% toolbox. BTL July 17, 2007.
%
% if called with an input string, that string will be appended to the message,
% eg IDfig('data from tf_120104_1.mat')
%
% if called with an output argument, it returns the handle to the text object it created
% 
% the tex interpreter is turned off, so that the _ characters don't result in subscripts
% BTL, Dec 12, 2004
%
% $Id: IDfig.m,v 1.3 2010/05/13 18:46:34 janosch Exp $
%
% modified by Rana   May 2010
ax = axis;
% figure out the name of the function which called this one
[st,ii] = dbstack;
if length(st)<2
    shortname = 'workspace';
else
    longname = st(2).name;
    index1 = max(find(longname == '\'));   % the last \ in the name
	if isempty(index1)
		shortname = longname;   % now compatible with v7 (doesn't return full path)
	else
		shortname = longname(index1+1:end);
    end
end
if nargin == 0
    bonus = [];
else
    bonus = [' ' varargin{1}];
end
% Add capability for a personal ident using AIRWOLF keyword
if strcmpi('airwolf',varargin{1})
  [osstat, osname] = system('uname'); % must use Linux or Mac OS
  if ~osstat
    if (strfind(osname,'Darwin') | strfind(osname,'Linux'))
      [mork, mindy] = system('whoami');    % get username
      [r2d2, c3p0]  = system('hostname');  % and hostname using system calls  
      % make the string, remembering to trim off the carraige return at the
      % end
      bonus = ['          by ' mindy(1:end-1) ' on ' c3p0(1:end-1)];
    end
  end
end
if strcmpi('kitt',varargin{1})
  [status,netstring] = system('net user');    % get username
  lines = textscan(netstring,'%s', 'delimiter', '\n');
  lines = strtrim(lines{:,:});
  username = char(lines(6));
  computer = char(lines(2));
  computer = strtrim(computer(21:end));
  % make the string, remembering to trim off the carraige return at the
  % end
  bonus = ['          by ' username ' on ' computer];
end
thing = text(ax(2), ax(3), ['created using ', shortname, '.m on ', date, bonus],...
    'VerticalAlignment','top','Rotation',90,'FontSize',10,'Interpreter','none',...
    'Color',[0.42 0.0 0.8]);
if nargout == 1
    tt = thing;
end
