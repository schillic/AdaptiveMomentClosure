function [approxBarMu,gamma]=derMatchClosure(net,mdyn)
% Auxiliary m-script used by momentClosure.m to compute the moment closure
% functions reported in the following reference.
%
% [1] Abhyudai Singh, João Hespanha. 
% LogNormal Moment Closures for Biochemical Reactions. 
% In Proc. of the 45th Conf. on Decision and Contr., Dec. 2006. 
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
% Created October 13, 2006 
%
% Modified June 27, 2007
% Get moment dynamics as a single structure

nMoments=size(mdyn.Mu.ndx,1);
nBarMu=size(mdyn.barMu.ndx,1);

C=zeros(nMoments,nMoments);
for s=1:nMoments
  for p=1:nMoments
    C(s,p)=pcom(mdyn.Mu.ndx(p,:),mdyn.Mu.ndx(s,:));
  end
end  

%C
gamma=zeros(nBarMu,nMoments);

for i=1:nBarMu
  c=zeros(nMoments,1);
  for s=1:nMoments
    c(s)=pcom(mdyn.barMu.ndx(i,:),mdyn.Mu.ndx(s,:));
  end  
  gamma(i,:)=(C\c)';  
end  

approxBarMu=sym(zeros(size(mdyn.barMu.ndx,1),1));
barMu=sym(zeros(size(mdyn.barMu.ndx,1),1));

for i=1:size(mdyn.barMu.ndx,1)
  approxBarMu(i)=simplify(prod(mdyn.Mu.sym.^(gamma(i,:)')));
end

% approxBarMu

%%%%%%%%%%%%%%%%%%%%

function c=pcom(m1,m2)

c=1;
for i=1:length(m1)
  c=c*gamma(m1(i)+1)/gamma(m2(i)+1)/gamma(m1(i)-m2(i)+1);
end  
