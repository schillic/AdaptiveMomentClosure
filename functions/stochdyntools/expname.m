function [name,texrules]=expname(net,ndx,prefix,latexprefix,latexsuffix)
% [name,texrules]=expname(net,ndx,prefix,latexprefix)
%
% Takes a list of monomials and produces a list of symbolic variables
% that correspond to the expected values of the monomials. For a given
% moment 
%        E [ X1^{m1} X2^{m2} X3^{m3} ]
% two formats are possible for the symbolic variable
%    (1)    {prefix}{m1}{m2}{m3}
% or
%    (2)    {prefix}_X1{m1}_X2{m2}_X3{m3}
% and two latex formats are available
%    (1)    {latexprefix}{m1}{m2}{m3}
% or 
%    (2)    {latexprefix}X1^{m1}X2^{m2}X3^{m3}{latexsuffix}
%
% The local variable 'format' decides which format to use. If the
% argument latexsuffix is present this variable is forced to
% 2, inluess latexperfix = Nan in which case it is forced to 1.

%
% Inputs
% ------
%
% net: Internal structure describing the network, typically
%      read from a .net file using  net=readNet(filename)
%
% ndx: Matrix describing the moments of interest, 
%      with one row per moment and one column per species. 
%      In particular, the m-th moment is given by
%         E [ X1^{m1} X2^{m2} X3^{m3} ... ]
%      where [a1 a2 a3 ...] is the m-th row of ndx.
%
%      The order in which species appear in the columns of ndxMu
%      is the same order in which they appear in net.Reaction,
%      which matches the order in which they were declared in
%      the .net file.
%
% prefix: optional string specifying the prefix used to create the
%         symbolic variables. By default 'mu' is used in both formats.
%
% latexprefix: optional string specifying the prefix used to create the
%         latex variables. By default '\mu' is used in format (A)
%         and '\E[' is used in format (B).
%
% latexsuffix: optional string specifying the prefix used to create the
%         latex variables. By default ']' is used in format (B)
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
% June 12, 2007
% Added support latex strings
%
% June 28, 2007
% Input is net structure

format=2;

if nargin<3
  prefix='mu';  
end  
  
if nargin<4
  latexprefix='\mu';  
  if format==2  
    latexprefix='\E[';  
  end  
end  
  
if nargin<5
  latexsuffix=']';  
else
  if isnan(latexsuffix)  
    format=1;
  else
    format=2;    
  end    
end  
  
nSpecies=length(net.species);
name=[];
texrules=[];
for i=1:size(ndx,1)
  if ~any(ndx(i,:))
    thisname='1';
    thistexname='1';
  else
    switch format
    case 1, % {prefix}{m1}{m2}{m3}   &   {latexprefix}{m1}{m2}{m3}
      thisname=[prefix,sprintf('_%d',ndx(i,:))];
      thistexname =[latexprefix,sprintf('_%d',ndx(i,:))];
    case 2,%{prefix}{m1}{m2}{m3}&{latexprefix}X1^{m1}X2^{m2}X3^{m3}{latexsuffix}
      thisname=prefix;   
      thistexname=latexprefix;      
      for j=1:nSpecies
	if ndx(i,j)==1
	  thisname=sprintf('%s_%s',thisname,net.species(j).id);
	  thistexname=sprintf('%s%s ',thistexname,net.species(j).idTex);
	elseif ndx(i,j)>1
	  thisname=sprintf('%s_%s%d',thisname,net.species(j).id,ndx(i,j));
	  thistexname=sprintf('%s%s^%d ',thistexname,net.species(j).idTex,ndx(i,j));
	end    
      end	
      assignin('base',thisname,sym(thisname));
      thistexname=[thistexname(1:end-1),latexsuffix];
    end
  end    
  name=1*[name;sym(thisname)]; % 1* prevent matrix([[]]) in scalars

  thisname=mylatex(sym(thisname));

  % quote '{}' in symbolic variable names to make it compatible with regexps  
  thisname=regexprep(thisname,'([\{\}])','\\$1');
    
  texrules=[texrules;{thisname,thistexname}];
end

if ~isempty(texrules)
  % reverse order to avoid replacement of partial names
  [dummy,k]=sort(sum(char(texrules(:,1))~=' ',2),1,'descend');
  texrules=texrules(k,:);
end

