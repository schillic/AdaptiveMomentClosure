function approxBarMu=zeroVarianceClosure(net,mdyn)
% Auxiliary m-script used by momentClosure.m to compute a moment closure
% functions that assumes zero variance, i.e. E[X^n]=E[X]^n.
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
% Created October 27, 2006 


nBarMu=size(mdyn.barMu.ndx,1);

symMoments=expname(net,eye(size(mdyn.barMu.ndx,2)));

approxBarMu=sym(zeros(size(mdyn.barMu.ndx,1),1));

%symMoments
%mdyn.barMu.ndx
for i=1:size(mdyn.barMu.ndx,1)
  approxBarMu(i)=simplify(prod(symMoments.^(mdyn.barMu.ndx(i,:)')));
end

%mdyn.barMu.ndx
%approxBarMu

