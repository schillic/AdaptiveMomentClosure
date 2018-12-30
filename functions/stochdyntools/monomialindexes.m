function [coef,mon,idx]=monomialindexes(poly,vars,mpnb)
% [coef,mon,idx]=monomialindexes(poly,vars)
%
% Returns the coefficients and indices of the monomials in a given
% polynomial
%
% Input:
% poly - polynomial on variables 'vars' (symbolica variable)
% vars - array with variables (array of symbolic variables)
%
% Output:
% coef - array of coefficients (symbolic variable)
% mon  - array of monomials (symbolic variable)
% idx  - indices of the monomial, one integers >=0 per variable
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
% Modified May 19, 2009
% Moved from maple to MuPAD
% 
% Modified Jul 19, 2009
% Optimized for MuPAD

vars=vars(:).'; % make vars a row vector
charvars=['[',regexprep(char(vars),'^matrix\(|\)$|\[|\]',''),']'];

if nargin<3
    mpnb=symengine;
end

%evalin(mpnb,['l:=poly2list(',char(poly),',',charvars,');']);
%coef=evalin(mpnb,'map(l,_index,1)');
%idx=double(evalin(mpnb,'matrix(map(l,_index,2))'));
%coef=feval(mpnb,'coeff',poly,charvars)
cmd=['[coeff(',char(poly),',',charvars,')];'];
coef=evalin(mpnb,cmd);
mon=feval(mpnb,'monomials',poly,charvars);
mon=mon./coef;
if isempty(mon)
    coef=mon;
end

if nargout>=3
  % convert symbolic expression of monomials into vectors of indexes
  cidx=feval(mpnb,'polylib::support',['poly(',char(poly),',',charvars,')']);

  if isempty(cidx)
      idx=[];
  else
      cidx=char(cidx);
      cidx=strrep(cidx,'matrix([','');
      cidx=strrep(cidx,'])','');
      cidx=strrep(cidx,'],','];');
      idx=eval(cidx);
      if size(idx,2) ~= length(vars)
          idx=idx';
      end
  end
end

if nargin<3
    %delete(mpnb)  % should be used, but apparently create errors
end

