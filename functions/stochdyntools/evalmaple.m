function in=evalmaple(in,old,new)

%evalmaple   Symbolic substitution.
%   evalmaple(S,OLD,NEW) replaces OLD with NEW in the symbolic expression S.
%   OLD is a symbolic variable, a string representing a variable name, or
%   a string (quoted) expression. NEW is a symbolic or numeric variable
%   or expression.
%
%   OLD and NEW are vectors or arrays of the same size, each element
%   of OLD is replaced by the corresponding element of NEW.  
%
% Attention: Unlike the evalsubs function, substitutions are not
% performed on already subtituted portions of the expression
%
% Copyright (C) 2006  Joao Hespanha

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% Modification history:
%
% Created October 28, 2008


if ~isempty(old)
  str='{';

  if ischar(old)
    if isnumeric(new)
      str=[old,'=',num2str(new)];
    else
      str=[old,'=',char(new)];
    end
  else
    for i=1:prod(size(old))
      if iscell(old)
	o=old{i};
      else
	o=old(i);
      end
      if iscell(new)
	n=new{i};
      else
	n=new(i);
      end
      if isnumeric(n)
	str=[str,char(o),'=',num2str(n),','];
      else
	str=[str,char(o),'=',char(n),','];
      end
    end
    str=str(1:end-1);
  end
  str=[str,'}'];
  
  if isa(in,'cell')
    for j=1:numel(in)
      if ischar(in{j})
	in{j}=maple('eval',in{j},str);
      else
	for i=1:numel(in{j})
	  in{j}(i)=maple('eval',char(in{j}(i)),str);
	end
      end
    end
  else
    if ischar(in)
      in=maple('eval',in,str);
    else
      for i=1:numel(in)
	in(i)=maple('eval',char(in(i)),str);
      end
    end
  end
end

