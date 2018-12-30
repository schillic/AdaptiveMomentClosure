function expr=subsParameters(net,expr,mpnb)
% new_expression=subsParameters(net,old_expression)
%
% This function takes a given expression and replaces the parameters
% from a network of chemical reactions by their default values.
%
% Inputs
% ------
%
% net    : Internal structure describing the network, typically
%          read from a .net file using  net=readNet(filename)
%
% old_expression : Vector of symbolic expressions where the
%          parameters should be eliminated
%
% Output
% ------
%
% new_expression : Vector of symbolic expressions obtained from
%          old_expression by replacing any parameters  by their
%          default values.
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
% Created November 4, 2006

if nargin<3
    mpnb=symengine;
end

old=[];
new=[];
for thisParameter=1:length(net.parameter)
  old=[old;cellstr(net.parameter(thisParameter).id)];
  new=[new;cellstr(char(sym(net.parameter(thisParameter).value)))];
end
expr=mysubs(expr,old,new,1,mpnb);

if nargin<3
    %delete(mpnb)  % should be used, but apparently create errors
end
