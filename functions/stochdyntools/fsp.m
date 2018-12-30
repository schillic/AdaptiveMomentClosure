function [newnet,newexpr]=fsp(oldnet,species,values,idTex,oldexpr)

mpnb=symengine;

verbose=0;

fprintf('converting ''%s'' to boolean... ',char(species));
t0=clock;

nSpecies=length(oldnet.species);
nReactions=length(oldnet.reaction);
nDerivatives=length(oldnet.raterule);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update unmodified sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newnet.parameter=oldnet.parameter;
newnet.substitutionrule=oldnet.substitutionrule;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update species section & create expression for # molecules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

found=0;
for thisSpecies=1:nSpecies
    if strcmp(oldnet.species(thisSpecies).id,char(species))
        if strcmp(oldnet.species(thisSpecies).type,'stochastic')
            found=1;
            break;
        else
            error('fsp: species ''%s'' is of type ''%s'', instead of ''stochastic''',species,net.species(thisSpecies).type);
        end
    end  
end

if ~found
    error('fsp: species ''%s'' not found',species);
end

newnet.species=[oldnet.species(1:thisSpecies-1);oldnet.species(thisSpecies+1:end)];

% create boolean species for given values
foundInitial=0;
speciesSym=sym(0);
for i=1:length(values)
    newnet.species(end+1,1).id=sprintf('%s_%d',oldnet.species(thisSpecies).id,values(i));
    speciesSym=speciesSym+values(i)*sym(newnet.species(end).id);
    newnet.species(end).idTex=sprintf('%s_%d',idTex,values(i));
    newnet.species(end).type='boolean';
    if values(i)==oldnet.species(thisSpecies).initialAmount
        newnet.species(end).initialAmount=1;
        foundInitial=1;
    else
        newnet.species(end).initialAmount=0;
    end
end

% create boolean species for remaining values
newnet.species(end+1,1).id=sprintf('%s_%d',oldnet.species(thisSpecies).id,NaN);
newnet.species(end).idTex=sprintf('%s_%d',idTex,NaN);
newnet.species(end).type='boolean';
newnet.species(end).initialAmount=1-foundInitial;

if verbose
    fprintf('old species\n');
    for i=1:length(oldnet.species)
        fprintf('\t%s\n',oldnet.species(i).id)
    end
    fprintf('new species\n');
    for i=1:length(newnet.species)
        fprintf('\t%s\n',newnet.species(i).id)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update derivatives section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newnet.raterule=oldnet.raterule;
for i=1:nDerivatives
    if strcmp(newnet.raterule(i).variable,char(species))
        error('fsp: derivatice not allowed for species ''%s''',char(species));
    end
    newnet.raterule(i).formula=char(mysubs(newnet.raterule(i).formula,sym(species),speciesSym,0,mpnb));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update reactions section
%%%%%%%%%%%%%%%%%%%%%%%%%%%

newnet.reaction=[];
for thisReaction=1:nReactions
    if verbose
        fprintf('processing reaction %d: rate = %s \n',thisReaction,oldnet.reaction(thisReaction).intensity);
        disp([oldnet.reaction(thisReaction).old,oldnet.reaction(thisReaction).new]);
    end

    % check is reaction needs to be expanded
    expand=0;
    vars=findsym(sym(oldnet.reaction(thisReaction).intensity));
    if strfind(vars,char(species))
        expand=1;
    else
        for j=1:length(oldnet.reaction(thisReaction).old)
            vars=findsym(sym(oldnet.reaction(thisReaction).old{j}));
            if strfind(vars,char(species))
                expand=1;
                break;
            end
            vars=findsym(sym(oldnet.reaction(thisReaction).new{j}));
            if strfind(vars,char(species))
                expand=1;
                break;
            end
        end
    end
    if expand==0
        % no need to expand
        newnet.reaction(end+1,1).intensity=oldnet.reaction(thisReaction).intensity;
        newnet.reaction(end,1).old=oldnet.reaction(thisReaction).old;
        newnet.reaction(end,1).new=oldnet.reaction(thisReaction).new;
        if verbose
            fprintf('  not expanding... new rate(%d) = %s \n',length(newnet.reaction),newnet.reaction(end).intensity);
            disp([newnet.reaction(end).old,newnet.reaction(end).new]);
        end
    else    
        % reaction needs to be expanded
        if verbose
            fprintf('  expanding...\n');
        end
        for i=1:length(values)
            intensity=char(sym(newnet.species(nSpecies+i-1).id)*mysubs(oldnet.reaction(thisReaction).intensity,sym(species),values(i),0,mpnb));
            if strcmp(intensity,'0')
                if verbose
                    fprintf('  new rate(skipping) = 0 (%s=%d)\n',char(species),values(i));
                end
                continue;
            end
            newnet.reaction(end+1,1).intensity=intensity;
            newnet.reaction(end).old=cell(0);
            newnet.reaction(end).new=cell(0);
            
            for j=1:length(oldnet.reaction(thisReaction).old)
                if strcmp(oldnet.reaction(thisReaction).old(j),species)
                    newval=mysubs(oldnet.reaction(thisReaction).new{j},sym(species),values(i),0,mpnb);
                    try
                        newval=double(newval);
                    catch
                        error('fsp: in reaction %d, new value of ''%s''=''%s'' is not a numerical value for a fixed value of this species',...
                              thisReaction,char(newval),char(species));
                    end
                    newnet.reaction(end).old{end+1,1}= ...
                        newnet.species(nSpecies+i-1).id;
                    newnet.reaction(end).new{end+1,1}='0';
                    k=find(values==newval);
                    if isempty(k)
                        k=length(values)+1;  % NaN
                    end
                    newnet.reaction(end).old{end+1,1}= ...
                        newnet.species(nSpecies+k-1).id;
                    newnet.reaction(end).new{end+1,1}='1';
                else
                    newnet.reaction(end).old{end+1,1}=oldnet.reaction(thisReaction).old{j};
                    newnet.reaction(end).new{end+1,1}=char(mysubs(oldnet.reaction(thisReaction).new{j},sym(species),values(i),0,mpnb));
                end 
           end
           if verbose
               fprintf('  new rate(%d) = %s  (%s=%d)\n',length(newnet.reaction),newnet.reaction(end).intensity,char(species),values(i));
               disp([newnet.reaction(end).old,newnet.reaction(end).new]);
           end
        end 
    end
end

newexpr=mysubs(oldexpr,sym(species),speciesSym,0,mpnb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add subsitution rules for incompatible boolean variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(values)+1
    for j=1:i-1
        if i~=j
            newnet.substitutionrule(end+1,1).old=[newnet.species(nSpecies+j-1).id,'*',newnet.species(nSpecies+i-1).id];
            newnet.substitutionrule(end,1).new='0';
        end
    end
end
if verbose
    for i=1:length(newnet.substitutionrule);
        fprintf('  rule(%d): {%s} > {%s}\n',i, ...
                newnet.substitutionrule(i).old,newnet.substitutionrule(i).new);
    end
end

        
fprintf('finished %.2fsec\n',etime(clock,t0));
