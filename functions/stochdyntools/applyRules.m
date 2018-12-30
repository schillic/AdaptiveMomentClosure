function in=applyRules(net,in,mpnb)
% out=applyRules(net,in)
% Applies the rules in net.substitutionrules to a symbolic expression
%
% Inputs
% ------
%
% net    : Internal structure describing the network, typically
%          read from a .net file using  net=readNet(filename)
%          All computations are done symbolic so any default values
%          stored in net.parameter are ignored.
% 
% in     : expression to be transformed
% 
% Output:
% -------
%
% out    : transformed expression
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
% Created May 8, 2007
%
% Modified June 26, 2007
% Updated structure field to net.substitutionrule instead of net.rule
%
% Modified Dec 28, 2008
% Using maple's 'applyrule' for recursive substitutions
%
% Modified May 19, 2009
% Moved from maple to MuPAD

if nargin<3
    mpnb=symengine;
end

nRules=length(net.substitutionrule);

if nRules>0
    [old{1:nRules,1}]=deal(net.substitutionrule(:).old);
    [new{1:nRules,1}]=deal(net.substitutionrule(:).new);

    in=expand(mysubs(expand(sym(in)),old,new,1,mpnb));
end

if nargin<3
    %delete(mpnb)  % should be used, but apparently create errors
end
