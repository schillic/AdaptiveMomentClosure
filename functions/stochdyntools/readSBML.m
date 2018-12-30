function [newSBMLModel,x0]=readSBML(filename)
% [net,x0]=readSBML(filename)
%
% This function reads an SBML file describing a network of chemical reactions
%
% Input
% ------
%
% filename : Name of the SBML file describing the chemical reactions
%
% Output
% ------
%
% net : Internal structure describing the network of chemical
%       reactions. This structure has two main fields: 
%         net.species    - describes the species involved 
%                          (from the 'species' section of the .net file)
%         net.parameters - describes the parameters that appear in
%                          the chemical reactions
%                          (from the 'parameters' section of the .net file)
%         net.reaction   - describes the chemical reactions
%                          (from the 'reactions' section of the .net file)
%
%       This retains all fields of the SBML structure returned by
%       TranslateSBML with the following data added
%
%         1) type of species (always stochastic)
%            net.species(i).type='stochastic';
%
%         2) parameters in net.reaction(j).kineticLaw added to net.parameter
%
%         3) fields added to describe reaction
%            net.reaction(i).old
%            net.reaction(i).new
%            net.reaction(i).intensity
%
% X0 : values for the SBML structure in species.initialAmounts 
%
% This is an incomplete function. Do not believe its output.
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

fprintf('WARNING (readSBML): This is an incomplete function. Do not believe its output.\n\t This function ignores all SBML units and assumes\n\t intensities in reaction per second\n')

SBMLModel=TranslateSBML(filename);
newSBMLModel=SBMLModel;

%% functionDefinitions not% supported
if ~isempty(SBMLModel.functionDefinition)
  error('readSBML: ''functionDefinition'' field not supported\n');
end

%% unitDefinitions not% supported
if ~isempty(SBMLModel.unitDefinition)
  warning('readSBML: ''unitDefinition'' field not supported\n');
end

%% only single compartment currently supported
if length(SBMLModel.compartment)~=1
  error('readSBML: only single compartment supported (%d compartments found)\n',length(SBMLModel.compartment));
end

%%%%%%%%%%
%% species
%%%%%%%%%%

x0=zeros(length(SBMLModel.species),1);
for thisSpecies=1:length(SBMLModel.species)
  % not SBML
  newSBMLModel.species(thisSpecies).type='stochastic';
  
  if SBMLModel.species(thisSpecies).boundaryCondition
    error('readSBML: species ''%s'' has boundary attribute, not yet implemented\n',SBMLModel.species(thisSpecies).id)
  end

  if SBMLModel.species(thisSpecies).constant
    error('readSBML: species ''%s'' has constant attribute, not yet implemented\n',SBMLModel.species(thisSpecies).id)
  end

  x0(thisSpecies)=SBMLModel.species(thisSpecies).initialAmount;
end

%%%%%%%%%%%%%
%% parameters
%%%%%%%%%%%%%

for thisParameter=1:length(SBMLModel.parameter)
  if ~SBMLModel.parameter(thisParameter).constant
    err=1;
    for thisRule=1:length(SBMLModel.rule)
      if strcmp(SBMLModel.rule(thisRule).typecode,'SBML_ASSIGNMENT_RULE') ...
	    & strcmp(SBMLModel.parameter(thisParameter).id, ...
		     SBMLModel.rule(thisRule).variable)
	err=0;
	break;
      end
    end
    if err
      fprintf('readSBML: parameter ''%s'' not constant, not yet implemented\n',...
	      SBMLModel.parameter(thisParameter).id)
    end
  end
end
for thisReaction=1:length(SBMLModel.reaction)
  for thisLaw=1:length(SBMLModel.reaction(thisReaction).kineticLaw)
    for thisParameter=1:length(SBMLModel.reaction(thisReaction).kineticLaw(thisLaw).parameter)
      newSBMLModel.parameter(end+1)=SBMLModel.reaction(thisReaction).kineticLaw(thisLaw).parameter(thisParameter);
      if ~SBMLModel.reaction(thisReaction).kineticLaw(thisLaw).parameter(thisParameter).constant
	err=1;
	for thisRule=1:length(SBMLModel.rule)
	  if strcmp(SBMLModel.rule(thisRule).typecode,'SBML_ASSIGNMENT_RULE') ...
		& strcmp(SBMLModel.parameter(thisParameter).id, ...
			 SBMLModel.rule(thisRule).variable)
	    err=0;
	    break;
	  end
	end
	if err
	  fprintf('readSBML: kineticLaw parameter ''%s'' not constant, not yet implemented\n',...
		  SBMLModel.reaction(thisReaction).kineticLaw(thisLaw).parameter(thisParameter).id)
	end
      end
    end
  end
end

%%%%%%%%%%%%
%% reactions
%%%%%%%%%%%%

for thisReaction=1:length(SBMLModel.reaction)
  
  old=[];
  new=[];
  
  %% stoichiometry
  
  % products
  for thisSpecies=1:length(SBMLModel.reaction(thisReaction).product)
    this=SBMLModel.reaction(thisReaction).product(thisSpecies).species;
    
    % check if species already in old
    k=1;
    while k<=length(old)
      if strcmp(old(k),this)
	break
      else
	k=k+1;
      end
    end
    if k>length(old)
      old=[old;cellstr(this)];
      new=[new;cellstr(this)];
    end
    new(k)=cellstr(char(simplify(sym(char(new(k)))-SBMLModel.reaction(thisReaction).product(thisSpecies).stoichiometry)));
    
    if SBMLModel.reaction(thisReaction).product(thisSpecies).denominator~=1
      fprintf('readSBML: reaction ''%s'', product ''%s'', denominator = %g not supported',...
	      SBMLModel.reaction(thisReaction).id,this,...
	      SBMLModel.reaction(thisReaction).product(thisSpecies).denominator);
    end
  end
  
  % reactants
  for thisSpecies=1:length(SBMLModel.reaction(thisReaction).reactant)
    this=SBMLModel.reaction(thisReaction).reactant(thisSpecies).species;
    
    % check if species already in old
    k=1;
    while k<=length(old)
      if strcmp(old(k),this)
	break
      else
	k=k+1;
      end
    end
    if k>length(old)
      old=[old;cellstr(this)];
      new=[new;cellstr(this)];
    end
    new(k)=cellstr(char(simplify(sym(char(new(k)))+SBMLModel.reaction(thisReaction).reactant(thisSpecies).stoichiometry)));
    
    if SBMLModel.reaction(thisReaction).reactant(thisSpecies).denominator~=1
      fprintf('readSBML: reaction ''%s'', reactant ''%s'', denominator = %g not supported',...
	      SBMLModel.reaction(thisReaction).id,this,...
	      SBMLModel.reaction(thisReaction).reactant(thisSpecies).denominator);
    end
  end
  newSBMLModel.reaction(thisReaction).old=old;
  newSBMLModel.reaction(thisReaction).new=new;
  
  %% rate
  
  if length(SBMLModel.reaction(thisReaction).kineticLaw)~=1
    error('readSBML: reaction ''%s'' has $d kinetic laws, not implemented\n',SBMLModel.reaction(thisReaction).id,length(SBMLModel.reaction(thisReaction).kineticLaw))
  end
  
  newSBMLModel.reaction(thisReaction).intensity= ...
      fixSyntax(SBMLModel.reaction(thisReaction).kineticLaw.formula); 
  
  try
    sym(newSBMLModel.reaction(thisReaction).intensity);
  catch
    error('readSBML: syntax error in reaction ''%s'' intensity ''%s''\n',SBMLModel.reaction(thisReaction).id,newSBMLModel.reaction(thisReaction).intensity) 
  end
  
  %% misc
  
  if SBMLModel.reaction(thisReaction).reversible
    error('readSBML: reaction ''%s'' has reversible attribute, not implemented\n',SBMLModel.reaction(thisReaction).id)
  end
  
  if SBMLModel.reaction(thisReaction).fast
    warning('readSBML: reaction ''%s'' has fast attribute, not yet implemented\n',SBMLModel.reaction(thisReaction).id)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rules - replace in reaction rates and stoichiometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for thisRule=1:length(SBMLModel.rule)
  if ~strcmp(SBMLModel.rule(thisRule).typecode,'SBML_ASSIGNMENT_RULE')
    fprintf('readSBML: %s = %s [%s], not yet implemented\n',...
	    SBMLModel.rule(thisRule).variable,fixSyntax(SBMLModel.rule(thisRule).formula),SBMLModel.rule(thisRule).typecode);
  end
  
  for thisReaction=1:length(SBMLModel.reaction)
    newSBMLModel.reaction(thisReaction).intensity=subsmaple(newSBMLModel.reaction(thisReaction).intensity,SBMLModel.rule(thisRule).variable,fixSyntax(SBMLModel.rule(thisRule).formula));
    newSBMLModel.reaction(thisReaction).new=subsmaple(newSBMLModel.reaction(thisReaction).new,SBMLModel.rule(thisRule).variable,fixSyntax(SBMLModel.rule(thisRule).formula));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking for polynomials resets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

checkPolynomial(newSBMLModel);

function str=fixSyntax(str)

% Applies a few rules to fix/simplify expressions

%str

% +-9 => -9
str=regexprep(str,'+-(\d)','-$1');  
% *power(XXX,-1)* => /XXX*
str=regexprep(str,'\*power\(([^\(]+),-1\)','/($1)');
str=regexprep(str,'^power\(([^\(]+),-1\)','1/($1)');

%str
