function approxBarMu=zeroCumulantsClosure(net,mdyn)
% Auxiliary m-script used by momentClosure.m to compute the moment closure
% functions obtained by setting high-order cumulants equal to zero.
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
% Created May 12, 2007
%
% Modified June 27, 2007
% Get moment dynamics as a single structure
%
% Modified May 19, 2009
% Moved from maple to MuPAD
%
% Modified August 25, 2010
% Prevent error when barMu is empty.

mpnb=symengine;

verbose=0;

nMoments=size(mdyn.Mu.ndx,1);
nBarMu=size(mdyn.barMu.ndx,1);


if nBarMu==0
    approxBarMu=sym([]);
    return
end

nSpecies=size(mdyn.barMu.ndx,2);

nSpecies=length(net.species);
nReactions=length(net.reaction);

%% create symbolic species vector

symSpecies=[];
for thisSpecies=1:nSpecies
  symSpecies=[symSpecies,sym(net.species(thisSpecies).id)];
end

% Compute polynomial expansion of the moment generating function, up
% to the order of the largest moment in BarMu
maxorder=max(sum(mdyn.barMu.ndx,2));
syms l__x
I=sym('sqrt(-1)');
g=taylor(exp(I*l__x),l__x,maxorder+1);
l__x=0;
for i=1:nSpecies
   l__x=l__x+sym(sprintf('l__%d',i))*symSpecies(i);
 end
%l__x
g=sym(subs(g,'l__x',l__x,0));
[coef,mon,idx]=monomialindexes(g,symSpecies,mpnb);

% take expected value
g=0;
for thisMon=1:size(idx,1)
    if any(idx(thisMon,:))
        %coef(thisMon)*expname(net,idx(thisMon,:))
        g=g+coef(thisMon)*expname(net,idx(thisMon,:));
    else
        g=g+coef(thisMon);
    end    
end  
%g


% find cumulants that must be set to zero:
% . moments that appear in BarMu
% . moments of non-maximal order that appear in moment generating function
% but not in mdyn.Mu.ndx

k=find(sum(idx,2)<maxorder);
zerocum=unique([mdyn.barMu.ndx;idx(k,:)],'rows');
for thisMon=size(zerocum,1):-1:1
  % monomial in state?    
  k=find(all(mdyn.Mu.ndx==ones(nMoments,1)*zerocum(thisMon,:),2));
  if ~any(zerocum(thisMon,:)) || ~isempty(k)
    zerocum(thisMon,:)=[];
  end    
end
%zerocum

% compute cummulants that must be set to zero
equation='';
unknown='';
for thisCum=1:size(zerocum,1)
  cg=log(g);
  for i=1:nSpecies
    cg=diff(cg,sprintf('l__%d',i),zerocum(thisCum,i));
  end
  %cg
  %[sprintf('l__%d=0,',1:nSpecies-1),sprintf('l__%d=0',nSpecies)]
  cum=expand(feval(symengine,'subs',cg,[sprintf('l__%d=0,',1:nSpecies-1),sprintf('l__%d=0',nSpecies)]));
  equation=[equation,',',char(cum)];
  unknown=[unknown,',',char(expname(net, zerocum(thisCum,:)))];
  %  expname(net, zerocum(thisCum,:))
end  
if length(equation)>1
  equation(1)=[];
  unknown(1)=[];
end  

% find moments
%equation
%unknown
sol=solve(equation,unknown);
if isempty(sol)
    error('Unable to find zeroCumulant moment closure\n');
end

% extract moment to output

approxBarMu=sym(zeros(size(mdyn.barMu.ndx,1),1));

for thisMon=1:size(mdyn.barMu.ndx,1)
  if any(mdyn.barMu.ndx(thisMon,:))
    if verbose  
      fprintf('cumulantsClosure: closing %s ~ ',char(expname(net,mdyn.barMu.ndx(thisMon,:))))
    end
    
    if isstruct(sol)
      approxBarMu(thisMon)=expand(getfield(sol,char(expname(net,mdyn.barMu.ndx(thisMon,:)))));
    else % single unknown?
      approxBarMu(thisMon)=expand(sol);
    end  
    
    if verbose  
      fprintf('%s\n',char(approxBarMu(thisMon)))
    end    
  else
    approxBarMu(thisMon)=1;
  end    
end

%mdyn.barMu.ndx
%approxBarMu

%delete(mpnb)  % should be used, but apparently create errors
