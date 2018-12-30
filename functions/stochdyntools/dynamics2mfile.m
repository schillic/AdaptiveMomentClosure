function dynamics2mfile(net,mdyn,funname,symParameters,preamble)
% dynamics2mfile(net,mdyn,funname,symParameters,preamble)
%
% Writes approximate moment dynamics to an m-file. Used by closureDynamics().
%
% Copyright (C) 2007  Joao Hespanha

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
% Created October 31, 2007 
% 
% Feb 27, 2010
%
% Added preamble

if nargin<4
    symParameters=[];
end

if nargin<5
    preamble=''
end

% make preamble

% add parameters  
if ~isempty(net.parameter)
  commented=0;    
  for i=1:length(net.parameter)
    include=1;
    for j=1:length(symParameters)
      if strcmp(char(symParameters(j)),net.parameter(i).id)
	include=0;
	break;
      end
    end
    if include
      if ~commented
	preamble=[preamble,sprintf('%%%% default parameters\n')];
	commented=1;	  
      end
      preamble=[preamble,sprintf('  %s = %s;\n',net.parameter(i).id,char(sym(net.parameter(i).value)))];
    end
  end
  % add constant variables
  commented=0;    
  for i=1:length(net.species)
    if strcmp(net.species(i).type,'constant')
      if ~commented
	preamble=[preamble,sprintf('%%%% constant species\n')];
	commented=1;	  
      end	  
      preamble=[preamble,sprintf('  %s = %s;\n',net.species(i).id,char(sym(net.species(i).initialAmount)))];
    end	
  end      
end

%preamble

if isfield(mdyn,'dotPhi')
  sym2mfile(funname,[mdyn.approxDotMu.sym;mdyn.dotPhi],[mdyn.Mu.sym;mdyn.Phi.sym],symParameters,preamble);
else    
  sym2mfile(funname,mdyn.approxDotMu.sym,mdyn.Mu.sym,symParameters,preamble);
end  

