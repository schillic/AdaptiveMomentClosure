function net2stochKit(net,filename,x0,symParameters,valueParameters)
% net2stochKit(net,filename,x0,symParameters,valueParameters)
%
% This function creates a C++ "ProgramDefinition" file that can be
% used to create a StochKit executable to run Monte Carlo
% simulations of a network of chemical reactions. 
%
% Attention: all the'post_reaction_counts' in the .net file must
% correspond to increments/decrements. 
%
% Input
% -----
%
% net      : Internal structure describing the network, typically
%          read from a .net file using  net=readNet(filename)
%
% filename : filename for the C++ "ProgramDefinition" file to be created. 
%
% x0       : vector of initial populations for the species.
%
% symParameters (optional) : Symbolic vector of parameters that
%         appear in the rate expressions (e.g., rate constants).
%
% valueParameters (optional): vector of numerical values for the symbolic
%                   parameters in symParameters 
%
% Output
% ------
%
% No output returned, but creates a C++ "ProgramDefinition" file
% that can be used to create a StochKit executables. To learn how
% to use this file, please consult the StochKit user guide at
%    http://www.engineering.ucsb.edu/~cse/StochKit/index.html
%
% ATTENTION: any existing m-file with the given name will be overwritten.
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
% Created October 28, 2006 
% 
% May 7, 2005
% Parameters can also be used to defined stoichiometry
%
% Jan 1, 2008
% rates converted to C code with cccode()

if nargin<5
  symParameters=[];
  valueParameters=[];
end

if length(net.raterule)>0
  error('StochKit cannot simulate continuous rates\n');
end

nSpecies=length(net.species);
nReactions=length(net.reaction);

fid=fopen(sprintf('%s.cpp',filename),'w');

%% Headers
fprintf(fid,'#include "ProblemDefinition.h"\n#include <stdlib.h>\n#include <iostream>\n\n');

%% Initialization
fprintf(fid,'Vector Initialize()\n{\n   Vector x0(%d,0.0);\n',nSpecies);
for i=1:length(x0)
  fprintf(fid,'   x0(%d) = %d;\n',i-1,round(x0(i)));
end
fprintf(fid,'   return x0;\n}\n\n');

%% Stoichiometry
fprintf(fid,'Matrix Stoichiometry()\n{\n   Matrix nu(%d,%d,0.0);\n\n',nSpecies,nReactions);

outputParameters(fid,net,symParameters,valueParameters);

for thisReaction=1:nReactions
  for i=1:length(net.reaction(thisReaction).old)
    % find species
    k=1;
    while k<=nSpecies
      if strcmp(char(net.reaction(thisReaction).old(i)),net.species(k).id)
	break;
      end
      k=k+1;
    end
    if k>nSpecies
      error('Species ''%s'' that appears in reaction %d is not knwon',...
	    char(net.reaction(thisReaction).old(i)),thisReaction),
    end
    change=char(sym(char(net.reaction(thisReaction).new(i)))-sym(char(net.reaction(thisReaction).old(i))));
%    if ~isnumeric(change)
%      error('Cannot find stoichiometry for species ''%s'' in reaction %d',...
%	    char(net.reaction(thisReaction).old(i)),thisReaction);
%    end
    fprintf(fid,'   nu(%d,%d) = %s;\n',k-1,thisReaction-1,change);
  end
end

fprintf(fid,'   return nu;\n}\n\n');

%% Propensity function
fprintf(fid,'Vector Propensity(const Vector& x)\n{\n   Vector lambda(%d);\n\n',nReactions);

outputParameters(fid,net,symParameters,valueParameters);

fprintf(fid,'// Species (defines)\n');
for i=1:nSpecies
  fprintf(fid,'#define %s x(%d)\n',net.species(i).id,i-1);
end
% rates
for i=1:nReactions
  cc=ccode(sym(net.reaction(i).intensity)); 
  fprintf(fid,'   lambda(%d) = %s;\n',i-1,cc(11:end));
end
fprintf(fid,'// Species (undefines)\n');
for i=1:nSpecies
  fprintf(fid,'#undef %s\n',net.species(i).id);
end

fprintf(fid,'\n   return lambda;\n}\n');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%
%% fprintf parameters
%%%%%%%%%%%%%%%%%%%%%

function outputParameters(fid,net,symParameters,valueParameters)

% defines for propensity functions -- net.parameter

comment=0;
if ~isempty(net.parameter)
  for i=1:length(net.parameter)
    include=1;
    % redefinition of net.parameter? if so do not define net.parameter
    for j=1:length(symParameters)
      if strcmp(char(symParameters(j)),net.parameter(i).id)
	include=0;
	break;
      end
    end
    if include
      if ~comment
	fprintf(fid,'// Parameters -- net.parameters\n');
	comment=1;
      end
      fprintf(fid,'   const float %s=%.9e;\n',net.parameter(i).id,double(net.parameter(i).value));
    end
  end
end

% defines for propensity functions -- symParameters
if ~isempty(symParameters)
  fprintf(fid,'// Parameters -- symParameters\n');
  comment=1;
  for i=1:length(symParameters)
    if ~isnumeric(valueParameters(i))
      error('Value for the symbolic parameter ''%s'' is not numeric (''%s'')',char(symParameters(i)),char(sym(valueParameters(i))))
    end
    fprintf(fid,'   const float %s=%.9e;\n',char(symParameters(i)),valueParameters(i));
  end
end

if comment
  fprintf(fid,'// Parameters -- end\n\n');
end
