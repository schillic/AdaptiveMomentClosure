function rc=checkPolynomial(net)
% rc=checkPolynomial(net)
%
% This function checks if the reaction rates and the post-reaction
% resets are polynomial functions of the molecule counts. 
%
% The stochastic formulation of moment dynamics according to
% Gillespie's model should always result in polynomial functions
% and some moment closure techniques require this to be so.
%
% Input
% ------
%
% net    : Internal structure describing the network, typically
%          read from a .net file using  net=readNet(filename)
%          All computations are done symbolic so any default values
%          stored in net.parameter are ignored.
%
% Output
% ------
%
% rc      : boolean variable that is true (=1) when everything is indeed
%           polynomials and false (=0) otherwise
%
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
% Created December 3, 2006 
%
% Modified May 5, 2007
% Fixed upper limit in reactions loop.
%
% Modified May 19, 2009
% Moved from maple to MuPAD

mpnb=symengine;

verbose=0;

% find all species symbols
symSpecies=[];
for thisSpecies=1:length(net.species)
  symSpecies=[symSpecies,sym(net.species(thisSpecies).id)];
end
% find all parameters symbols
symParameters=[];
for thisParameters=1:length(net.parameter)
  symParameters=[symParameters,sym(net.parameter(thisParameters).id)];
end

rc = 1;
for i=1:length(net.reaction)
  if verbose || (~mod(i,10) && length(net.reaction)>20)
    fprintf('checkPolynomial: checking reaction %d, intensity: ''%s''...\n',i,net.reaction(i).intensity);
  end
  rc = rc & check1poly(symSpecies,symParameters,net.reaction(i).intensity,'intensity',mpnb);
  for j=1:length(net.reaction(i).old)
      if all(symSpecies~=sym(char(net.reaction(i).old(j))))
	rc = 0;
	error('checkPolynomial: ''%s'' in reaction %d list of species is not a valid species',net.reaction(i).old(j),i)
      end
  end
  for j=1:length(net.reaction(i).new)
    if verbose
      fprintf('checkPolynomial: checking reaction %d, post-reaction species %d: ''%s''\n',i,j,net.reaction(i).new{j});
    end
    rc = rc & check1poly(symSpecies,symParameters, ...
			 net.reaction(i).new(j),'post-reaction value',mpnb);
  end
end

if ~rc
  fprintf('\nATTENTION: some moment closure functions may not work since not all intesities and/or resets are polynomials on the populations, with coefficients depending on the parameters\n\n');
end

%delete(mpnb)  % should be used, but apparently creates errors


function rc=check1poly(symSpecies,symParameters,expr,str,mpnb)

rc=1;

try
  coef=monomialindexes(sym(char(expr)),symSpecies,mpnb);
catch
  fprintf('checkPolynomial: the following %s is not a polynomial on the species\n\t %s: ''%s''\n',str,str,char(expr));
  fprintf('\t species: ''%s''\n',findsym(symSpecies))
  rc=0;
  return 
end

coefsym=sym(['[',findsym(coef),']']);
for i=1:length(coefsym)
  if all(symParameters~=coefsym(i))
    fprintf('checkPolynomial: unknown coefficient ''%s'' in the polynomial ''%s''\n',char(coefsym(i)),char(expr));
    rc=0;
  end
end

