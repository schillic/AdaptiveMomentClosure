function net=readNet(filename)
% net=readNet(filename)
%
% This function reads a .net file describing a network of chemical reactions
%
% Input
% ------
%
% filename : Name of the file describing the chemical reactions
%
% Output
% ------
%
% net : Internal structure describing the network of chemical
%       reactions. This structure has four main fields: 
%         net.species    - describes the species involved 
%                          (from the 'species' section of the .net file)
%         net.parameters - describes the parameters that appear in
%                          the chemical reactions
%                          (from the 'parameters' section of the .net file)
%         net.substitutionrule - describes the algebraic rules that
%                          should be used to simplify computations
%                          that involve populations 
%                          (from the 'rules' section of the .net file)
%         net.reaction   - describes the chemical reactions
%                          (from the 'reactions' section of the .net file)
%         net.raterule   - describes continuous rates of change for
%                          the populations  
%                          (from the 'derivatives' section of the .net file)
%
% .net file format
% ----------------
%
% The file describing the chemical reactions has the following format
%
% The 'species' section describes the chemical species involved in
% the network. It also specifies assumptions that can be made in
% analysing the system.  
%
% Following the 'species:' header, each line is of the form:
%
% species_name "species_latex_name" species_type initial_value;
%
%   'species_name' stands for the symbol that identifies a chemical
%   species. This symbol should contain no spaces and be a valid name
%   for a Matlab variable.
%
%   'species_latex_name' is an optional a latex string that represents the
%   species population. This string appears quoted by " and should not
%   include any other quotes ("). 
%
%   'species_type' stands for the assumed type for the species. 
%
%   'initial_value' is an optional number that specifies the initial
%   number of molecules for numerical simulations. 
%
%   The 'species type' specifies assumptios that can be made in
%   analysing the system. It can be one of the following keywords:
%    'boolean'  - only 0 or 1 molecules can be present
%    'constant' - species population remains constant
%    'deterministic' - the expected number of molecules is much
%           larger than its standard deviation and therefore
%           stochastic effects can be negleted for this species
%    'stochastic' - no assumptios should be made about this species
%
% The (optional) 'parameters' section provides default numerical
% values for parameters that appear in the 'reactions' section.
%
% Following the 'parameters:' header, each line is of the form:
%
% parameter_name "parameter_latex_name" = default_value;
%
%   'parameter_name' stands for the symbol that identifies the
%   parameter. This symbol should contain no spaces and be a valid
%   name for a Matlab variable.
%
%   'parameter_latex_name' is an optional a latex string that represents the
%   parameter. This string appears quoted by " and should not include
%   any other quotes (").
%
%   'default_value' is a symbolic expression that stands for the default
%   value of the parameter. All symbolic computations ignore this value
%   and treat parameters as symbolic variables. Default values are only
%   used when numerical values are needed.
%
% The (optional) 'rules' section provides substitution rules that
% should be used to simplify computations.
%
% Following a 'rules:', each line is of the form: 
%
% {old_expression} > {new_expression};
%
%   'old_expression' stands for a valid Matlab expression that should
%   be replaced by 'new_expression' 
%
% The 'reactions' section describes the chemical reactions involved in
% the network, including their stochiometry and rates.
%
% Following the 'reactions:' header, each line is of the form:
% 
% rate = rate_expr;  {list_species} > {post_reaction_counts};
%
%   'rate_expr' stands for an expression describing the rate at
%   which the reaction occurs. It should be a polynomial expression.
%
%   'list_species' is a comma-separated list with the symbols of the
%   chemical species whose stochiometry changes in the reaction
%
%   'post_reaction_counts' is a comma-separated list of expressions
%   that specifies how the molecule counts for each symbol changes.
%
% The 'derivatives' section describes the continuous evolutions for the
% populations of the species.
%
% Following the 'derivatives:' header, each line is of the form: 
% 
% d/dt species += derivative_expr;
%
%   'species' is the symbol of the chemical species whose derivative is
%   provided 
%
%   'derivative_expr' is an expression describing the rate of change of
%   the population. This rate is to be 'added' to the discrete rules
%   specified in the 'reactions' section and to other derivatives for the
%   same species specified in the 'derivatives' section.
%
% Comments can be inserted anywhere with the prefix '%'
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
% November 4, 2006
% Added support for net.parameter
%
% June 12, 2007
% Added support latex strings
%
% June 27, 2007
% Added support for derivative rules

verbose=0;

t0=clock;
fprintf('reading file ''%s''... ',filename);

% read & parse file
[rawLines,statements,statementLines]=textparse(filename);
%statements

stat=1;

%%%%%%%%%%%%%%%%%%%%%%
% initialize structure
%%%%%%%%%%%%%%%%%%%%%%

net.species=[];
symSpecies=[];

net.parameter=[];
symParameters=[];

net.substitutionrule=[];

net.reaction=[];

net.raterule=[];

%%%%%%%%%%%%%%%%%
% reading species
%%%%%%%%%%%%%%%%%

if ~strcmp(lower(statements{stat}),'species')
  error('readNet: syntax error in line %d\n>>>%s\nreadNet: keyword ''species'' expected, instead found ''%s''',...
	statementLines(stat),rawLines{statementLines(stat)},char(statements{stat}))  
end  

while 1
  stat=stat+1;

  % 'species' format  
  % species_name ["species_latex_name"] species_type [initial_value];
  tokens=regexp(statements{stat},'^(\w+)\s+("[^"]*")?\s*(\w+)\s*(.*)$','tokens','once');

  if length(tokens)==4 & isstr(tokens{1}) & isstr(tokens{3})
    if ~checkId(tokens{1})
        error('readNet: invalid species name ''%s''\n',tokens{1});
    end
    net.species(end+1,1).id=tokens{1};
    if isempty(tokens{2})
      net.species(end).idTex=tokens{1};
    else      
      net.species(end).idTex=tokens{2}(2:end-1);
    end
    net.species(end).type=tokens{3};
    if isempty(tokens{4})
      tokens{4}='NaN';
      warning('\nreadNet: no initial value for  species ''%s'' Expect problems\n',tokens{1});
    end
  else
    break;
  end

  if ~any(ismember({'stochastic','boolean','deterministic','constant'},tokens{3}))
    error('readNet: unknown type ''%s'' for species ''%s''\n',tokens{3},tokens{1});
  end    
  
  try
    symSpecies=[symSpecies,sym(net.species(end).id)];
    net.species(end).initialAmount=sym(tokens{4});
  catch
    error('readNet: syntax error in line %d\n>>>%s\nreadNet: invalid species name ''%s''',...
	  statementLines(stat),rawLines{statementLines(stat)},net.species(end).id) 
  end
  
end

if isempty(net.species)
    error('readNet: syntax error in line %d\n>>>%s\nreadNet: species name and type expected, instead found ''%s''',...
	  statementLines(stat),rawLines{statementLines(stat)},char(statements{stat})) 
end


if verbose
  disp('net.species(:)')
  for i=1:length(net.species)
    disp(net.species(i))
    fprintf('    initialAmount= %s\n\n',char(net.species(i).initialAmount));
  end
end

%%%%%%%%%%%%%%%%%%%%
% reading parameters
%%%%%%%%%%%%%%%%%%%%

if strcmp(lower(statements{stat}),'parameters')
  while 1
    stat=stat+1;

    % 'parameters' format
    % parameter_name = default_value;
    tokens=regexp(statements{stat},'^(\w+)\s*("[^"]*")?\s*=\s*(.*)$','tokens','once');

    if length(tokens)==3
        if ~checkId(tokens{1})
            error('readNet: invalid parameter name ''%s''\n',tokens{1});
        end
      net.parameter(end+1,1).id=tokens{1};
      if isempty(tokens{2})
	net.parameter(end).idTex=tokens{1};
      else
	net.parameter(end).idTex=tokens{2}(2:end-1);
      end      
      net.parameter(end).value=tokens{3};
    else
      break;
    end
  
    try
      symParameters=[symParameters,sym(net.parameter(end).id)];
    catch
      error('readNet: syntax error in line %d\n>>>%s\nreadNet: invalid parameter name ''%s''',...
	  statementLines(stat),rawLines{statementLines(stat)},net.parameter(end).id) 
    end
    assignin('base',net.parameter(end).id,sym(net.parameter(end).id));

    try
      sym(net.parameter(end).value);
    catch
      error('#line %d: ''%s''\nreadNet: syntax error in default parameter value ''%s''',...
	  statementLines(stat),rawLines{statementLines(stat)},net.parameter(end).value);
    end
    net.parameter(end).value=sym(net.parameter(end).value);
    
  end
end

if verbose
  disp('net.parameter(:)')
  for i=1:length(net.parameter)
    disp(net.parameter(i))
    fprintf('    value= ''%s''\n\n',char(net.parameter(i).value));
  end
end

%%%%%%%%%%%%%%%%%%%%
% reading rules
%%%%%%%%%%%%%%%%%%%%

if strcmp(lower(statements{stat}),'rules')
  while 1
    stat=stat+1;
    
    % 'rules' format
    % {old_expression} > {new_expression};
    tokens=regexp(statements{stat},'^{(.*)}\s*>\s*{(.*)}$','tokens','once');

    if length(tokens)==2
        net.substitutionrule(end+1,1).old=tokens{1};
        net.substitutionrule(end).new=tokens{2};
    else
      break;
    end
  
    try
      sym(net.substitutionrule(end).old);
    catch
      error('readNet: syntax error in line %d\n>>>%s\nreadNet: invalid rule ''old'' field ''%s''',...
	  statementLines(stat),rawLines{statementLines(stat)},net.substitutionrule(end).old) 
    end
    
    try
      sym(net.substitutionrule(end).new);
    catch
      error('readNet: syntax error in line %d\n>>>%s\nreadNet: invalid rule ''old'' field ''%s''',...
	    statementLines(stat),rawLines{statementLines(stat)},net.substitutionrule(end).old)
    end
    
  end
end

if verbose
  disp('net.substitutionrule(:)')
  for i=1:length(net.substitutionrule)
    net.substitutionrule(i)
  end
end

%%%%%%%%%%%%%%%%%%%
% reading reactions
%%%%%%%%%%%%%%%%%%%

if ~strcmp(lower(statements{stat}),'reactions')
  error('#line %d: ''%s''\nreadNet: unexpected statement ''%s''.\n\t expecting keywords ''reactions''',...
	statementLines(stat),rawLines{statementLines(stat)},char(statements{stat}))  
end  

while stat<length(statements)
  stat=stat+1;

  % 'reactions' format
  % rate = rate_expr;  {list_species} > {post_reaction_counts};
  % parse: rate = rate_expr;
  tokens=regexp(statements{stat},'^rate\s*=\s*(.*)$','tokens','once');

  if length(tokens)==1
    net.reaction(end+1,1).intensity=tokens{1};
  else  
    break;    
  end    
  
  % check for syntax error
  try
    sym(net.reaction(end,1).intensity);
  catch
    error('#line %d: ''%s''\nreadNet: syntax error in rate ''%s''',...
	statementLines(stat),rawLines{statementLines(stat)},net.reaction(end,1).intensity);
  end
  
  stat=stat+1;
  if stat>length(statements)
      error('#line %d: ''%s''\nreadNet: unexpected end of file',...
	statementLines(stat-1),rawLines{statementLines(stat-1)});
  end

  % parse: {list_species} > {post_reaction_counts};
  tokens=regexp(statements{stat},'^{(.*)}\s*>\s*{(.*)}$','tokens','once');

  if length(tokens)==2
    net.reaction(end).old=list2char(tokens{1});
    net.reaction(end).new=list2char(tokens{2});
  else  
    error('#line %d: ''%s''\nreadNet: ''{list} = {list}'' expected, instead found ''%s''',...
	  statementLines(stat),rawLines{statementLines(stat)},char(statements{stat}))
  end

  % check for syntax error
  for ii=1:length(net.reaction(end).old)
    if all(symSpecies~=sym(char(net.reaction(end).old(ii))))
      error('#line %d: ''%s''\nreadNet: symbol ''%s'' in pre-reaction list is not a valid species\n\t species: ''%s''\n',...
	    statementLines(stat),rawLines{statementLines(stat)},char(net.reaction(end).old(ii)),findsym(symSpecies));
    end
  end
  
  try
    for ii=1:length(net.reaction(end).new)
      sym(char(net.reaction(end).new(ii)));
    end
  catch
    error('syntax error in post-reaction ''%s''',...
	  statementLines(stat),rawLines{statementLines(stat)},char(net.reaction(end).new(ii)));
  end
  
end  

if verbose
  disp('net.reaction(:)')
  for i=1:length(net.reaction)
    disp(net.reaction(i))
  end
end

%%%%%%%%%%%%%%%%%%%
% reading rates
%%%%%%%%%%%%%%%%%%%

if stat<length(statements) & ~strcmp(lower(statements{stat}),'derivatives')
  error('#line %d: ''%s''\nreadNet: unexpected statement ''%s''.\n\t expecting ''derivatives''',...
	statementLines(stat),rawLines{statementLines(stat)},char(statements{stat}))  
end  

while stat<length(statements);
  stat=stat+1;

  % 'rates' format
  % d/dt species += derivative_expr;
  tokens=regexp(statements{stat},'^d/dt\s*(\w+)\s*\+\=\s*(.*)$','tokens','once');

  if length(tokens)==2
    net.raterule(end+1,1).variable=tokens{1};
    net.raterule(end,1).formula=tokens{2};
  else  
    error('#line %d: ''%s''\nreadNet: ''d/dt species += expr'' expected, instead found ''%s''',...
	  statementLines(stat),rawLines{statementLines(stat)},char(statements{stat}))
  end    
  
  % check for syntax error
  try
    sym(net.raterule(end,1).variable);
  catch
    error('#line %d: ''%s''\nreadNet: syntax error in species ''%s''',...
	statementLines(stat),rawLines{statementLines(stat)},net.raterule(end,1).variable);
  end
  
  try
    sym(net.raterule(end,1).formula);
  catch
    error('#line %d: ''%s''\nreadNet: syntax error in derivative expression ''%s''',...
	statementLines(stat),rawLines{statementLines(stat)},net.raterule(end,1).formula);
  end
  
  if all(symSpecies~=sym(char(net.raterule(end).variable)))
    error('#line %d: ''%s''\nreadNet: symbol ''%s'' in species is not a valid species\n\t species: ''%s''\n',...
	statementLines(stat),rawLines{statementLines(stat)},char(net.raterule(end).variable),findsym(symSpecies));
  end
  
end  

if verbose
  disp('net.raterule(:)')
  for i=1:length(net.raterule)
    disp(net.raterule(i))
  end
end

fprintf('finished %.2fsec\n',etime(clock,t0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking for polynomials resets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

checkPolynomial(net);


function list=list2char(str)

k=find(str==' ');str(k)=[];  % remove spaces
list=[];
k=find(str==',');
s=1;
for i=k
  list=[list;cellstr(str(s:i-1))];
  s=i+1;  
end  
list=[list;cellstr(str(s:end))];

%%%% Check if species/parameter id is valid

function rc=checkId(str)

if strcmp(str,'E')
    fprintf('\nreadNet: species name ''%s'' conflicts with mupad protected symbol\n',str);
    rc=false;
else
    rc=true;
end
    
