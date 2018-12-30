function Lpsi=applyGenerator(net,psi,mpnb)
% Lpsi=applyGenerator(net,psi)
%
% Apply to the vector-valued function psi the infinitesimal
% generator of the Markov process defined by network of chemical
% reactions.
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
% Created March 29. 2010 from code previously in momentDynamics.m
    
verbose=0;

MuPADloops=1;  % do more computations inside MuPAD
if MuPADloops && nargin<3
    mpnb=symengine;
end

nReactions=length(net.reaction);
nDerivatives=length(net.raterule);
nMoments=length(psi);

if MuPADloops
    %%%% Create variables with intensities and reset maps in MuPAD
    intensity=sprintf('%s,',net.reaction(:).intensity);
    intensity=['intensity:=[',intensity(1:end-1),']'];
    
    clear old new
    [old{1:nReactions,1}]=deal(net.reaction(:).old);
    old=cellfun(@(x) (['[',sprintf('%s,',x{:}),']']),old,'UniformOutput',0);
    old=regexprep(['olds:=[',sprintf('%s,',old{:}),']'],',]',']');
    [new{1:nReactions,1}]=deal(net.reaction(:).new);
    new=cellfun(@(x) (['[',sprintf('%s,',x{:}),']']),new,'UniformOutput',0);
    new=regexprep(['news:=[',sprintf('%s,',new{:}),']'],',]',']');
    evalin(mpnb,intensity);
    evalin(mpnb,old);
    evalin(mpnb,new);
end


Lpsi=sym(zeros(nMoments,1));

for thisMoment=1:nMoments
    if verbose || ~mod(thisMoment,10)
        fprintf('\n   ...applyGenerator : computing moment %d/%d... working... ',thisMoment,nMoments);
    end
    
    %%% generator formula - loop over reactions
    if MuPADloops
        %%%% computation in MuPAD
        expression=char(psi(thisMoment));
        cmd=sprintf('der:=0;_for(k,1,nops(intensity),1,(der:=expand(der+intensity[k]*(subs(%s,zip(olds[k],news[k],_equal))-%s))));',expression,expression);
        Lpsi(thisMoment)=evalin(mpnb,cmd);
    else
        %%%% computations in MATLAB
        Lpsi(thisMoment)=0;
        for thisReaction=1:nReactions
            Lpsi(thisMoment)=expand(Lpsi(thisMoment)+sym(net.reaction(thisReaction).intensity)*(mysubs(psi(thisMoment),net.reaction(thisReaction).old,net.reaction(thisReaction).new,0,mpnb)-psi(thisMoment)));
            if 0
                psi(thisMoment)
                net.reaction(thisReaction).intensity
                [net.reaction(thisReaction).old,net.reaction(thisReaction).new]
                
                [psi(thisMoment),mysubs(psi(thisMoment),net.reaction(thisReaction).old,net.reaction(thisReaction).new,0,mpnb)]
                
                expand(Lpsi(thisMoment))
            end
        end  % for thisReaction
    end
    
    %%% generator formula - loop over derivatives
    for thisDerivative=1:nDerivatives
        Lpsi(thisMoment)=expand(Lpsi(thisMoment)+diff(psi(thisMoment),net.raterule(thisDerivative).variable)*net.raterule(thisDerivative).formula);
        
        if 0
            psi(thisMoment)
            net.raterule(thisDerivative).variable
            net.raterule(thisDerivative).formula
            
            diff(psi(thisMoment),net.raterule(thisDerivative).variable)
            
            Lpsi(thisMoment)
        end
    end  % for thisDerivative
    
end % for thisMoment

%% apply substitution rules
Lpsi=applyRules(net,Lpsi);

if MuPADloops && nargin<3
    %    delete(mpnb)
end

