function [muEq,phiEq,xiEq]=momentEquilibrium(net,mdyn,numerical,nophi,smallcv)
% [muEq,phiEq,xiEq]=momentEquilibrium(net,mdyn,numerical,nophi,smallcv)
%
% This function computes the equilibrium points of given
% approximate (closed) moment dynamics, of the form 
%           \dot Mu = A Mu + B F(\Mu)      ($$)
%     where 
%       A,B   - appropriate matrices from exact moment dynamics
%       Mu    - column vector that approximates Mu
%       F()   - moment closure function
% or 
%           \dot Phi = G(Phi)
%           \dot Mu  = A Mu + B F(Phi,Mu)  ($$$)
%     where 
%       Mu    - column vector that approximates Mu
%       F()   - moment closure function
%       Phi   - state of an additional dynamical system
%
% Inputs
% ------
%
% net is a structure describing the network of chemical reactions.
%     It is typically obtained from a .net file using net=readNet(filename) 
%
% mdyn is a structure characterizing the moment dynamics, typically
%      computed using closureDynamics(), momentDynamics(), or momentClosure(). 
%      Details can be obtained with 'help mdyn'
%
% numerical (optional) when nonzero the parameters should be replaced
%           by their numerical values prior to finding the
%           equilibrium. If 'numerical' is negative then only real, 
%           positive solutions are returned. 
%
% nophi (optional) when nonzero the function does NOT replace the
%           equilibrium point for the Phi subsystem in ($$$) into the
%           equilibrium point for the Mu system.
%
% smallcv (optional) when nonzero an approximate equilibrium is
%         computed assuming that the coefficients of variation are
%         much smaller than one for all variables. Currently does not
%         work for Van Kampen's linear noise approximation 
%
% Outputs
% -------
%
% muEq is cell array of symbolic vectors with the equilibrium points
%      of the vector Mu. Each entry of the cell array corresponds to 
%      one equilibrium point.
%
% phiEq is a cell array of symbolic vectors with the equilibrium
%       points of first-order "perturbation" vector Xi in for the
%       linear noise approximation (see PDF documentation). Each entry of 
%       the cell array corresponds to one equilibrium point. 
%       This output is only computed for approximate dynamics of the
%       form ($$$) 
%
% xiEq is a cell array of symbolic vectors with the equilibrium points
%      of the perturbations Xi around Phi.Each entry of the cell array
%      corresponds to one equilibrium point. 
%      This output is only computed for the linear noise approximation.
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
% Created December 25, 2007 
%
% January 10, 2007
% updated inline documentation

if nargin <3
  numerical=0;
end

if nargin <4
  nophi=0;
end

if nargin <5
  smallcv=0;
end

mpnb=symengine;

phiEq=[];
xiEq=[];

fprintf('computing moment equilibrium (numerical=%d,nophi=%d,smallcv=%d)... ',numerical,nophi,smallcv);
t0=clock;

if numerical
  mdyn.approxDotMu.sym=subsParameters(net,mdyn.approxDotMu.sym,mpnb); 
  if isfield(mdyn,'dotPhi')
    mdyn.dotPhi=subsParameters(net,mdyn.dotPhi,mpnb);
  end    
  if isfield(mdyn,'approxDotXi')
    mdyn.approxDotXi=subsParameters(net,mdyn.approxDotXi,mpnb);
  end    
  if isfield(mdyn.Mu,'xi')
    mdyn.Mu.xi=subsParameters(net,mdyn.Mu.xi,mpnb);
  end    
end

if ~smallcv

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% "Exact" computation (no smallcv)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% solve deterministic equations
    if isfield(mdyn,'dotPhi') && (nargin>=2 || ~nophi)
        phi_equations=sym2cell(mdyn.dotPhi);
        phi_variables=sym2cell(mdyn.Phi.sym);
        phi_solution=solve(phi_equations{:},phi_variables{:});
        % typically multiple solutions  
        phiEq=subssolve(mdyn.Phi.sym,phi_variables,phi_solution,mpnb);
    end
    
    %% solve for perturbations
    if isfield(mdyn,'approxDotXi') && nargin>=3
        xi_equations=sym2cell(mdyn.approxDotXi.sym);
        xi_variables=sym2cell(mdyn.Xi.sym);
        xi_solution=solve(xi_equations{:},xi_variables{:});
        xiEq=subssolve(mdyn.Xi.sym,xi_variables,xi_solution,mpnb);
        % single solution  
        xiEq=xiEq{1}; 
    end
    
    %% solve Mu equilibrium
    equations=sym2cell(mdyn.approxDotMu.sym);
    variables=sym2cell(mdyn.Mu.sym);
    solution=solve(equations{:},variables{:});
    muEq=subssolve(mdyn.Mu.sym,variables,solution,mpnb);
    
    if isfield(mdyn,'dotPhi') && ~nophi
        % single solution
        muEq=subssolve(muEq{1},phi_variables,phi_solution,mpnb);
    end

else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Abhi's linear noise approximation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% find 1st and 2nd order moments
    
    order=sum(mdyn.Mu.ndx,2);
    nSpecies=length(net.species);
    first=find(order==1);
    
    if (any(any(mdyn.Mu.ndx(first,:)~=eye(nSpecies))))
        error('momentDynamics is producing incompatible Mu vectors')
    end  
    
    %% Create nominal values
    [mdyn.Phi.sym,texrules1]=expname(net,mdyn.Mu.ndx(first,:),'phi','\phi_{','}');
    
    %% Create perturbations Xi
    nMoments=size(mdyn.Mu.ndx,1);
    [mdyn.Xi.sym,texrules2]=expname(net,mdyn.Mu.ndx,'xi','\xi_{','}');
    texrules=[texrules2;texrules1];
    
    %mdyn.Mu.sym
    %mdyn.Phi.sym
    %mdyn.Xi.sym
    
    mdyn.Xi.mu=sym(zeros(nMoments,1));
    % first order moments
    for i=first'
        mdyn.Xi.mu(i)=mdyn.Mu.sym(i)/mdyn.Phi.sym(i)-1;
    end
    % higher order moments
    for i=find(order>1)'
        mdyn.Xi.mu(i)=mdyn.Mu.sym(i)/prod(mdyn.Mu.sym(first).^(mdyn.Mu.ndx(i,:)'))-1;
    end
    
    %mdyn.Xi.mu  
    
    %% Express dynamics in terms of perturbations Xi
    
    equations=sym2cell(mdyn.Xi.mu-mdyn.Xi.sym);
    variables=sym2cell(mdyn.Mu.sym);
    solution=solve(equations{:},variables{:});
    
    % xx=struct2cell(solution);xx{:}
    
    muEq=subssolve(mdyn.approxDotMu.sym,variables,solution,mpnb);
    mdyn.Mu.xi=subssolve(mdyn.Mu.sym,variables,solution,mpnb);
    mdyn.Mu.xi=mdyn.Mu.xi{1};
    %  mdyn.Mu.xi
    approxDotMu=expand(muEq{1});
    
    %% linear expansion
    
    approxDotMu0 = mysubs(approxDotMu,mdyn.Xi.sym,zeros(nMoments,1),0,mpnb);
    approxDotMu1 = mysubs(jacobian(approxDotMu,mdyn.Xi.sym),mdyn.Xi.sym,zeros(nMoments,1),0,mpnb)*mdyn.Xi.sym;
    
    %% solve deterministic equations
    
    phi_equations=sym2cell(approxDotMu0(first));
    phi_variables=sym2cell(mdyn.Phi.sym);
    phi_solution=solve(phi_equations{:},phi_variables{:});
    
    %% solve for perturbations
    
    xi_equations=sym2cell(approxDotMu0+approxDotMu1);
    xi_variables=sym2cell(mdyn.Xi.sym);
    xi_solution=solve(xi_equations{:},xi_variables{:});
    
    %% compute equilibrium
    
    % typically multiple solutions  
    phiEq=subssolve(mdyn.Phi.sym,phi_variables,phi_solution,mpnb);
    
    xiEq=subssolve(mdyn.Xi.sym,xi_variables,xi_solution,mpnb);
    % single solution
    xiEq=xiEq{1}; 
    
    muEq=subssolve(mdyn.Mu.xi,xi_variables,xi_solution,mpnb);
    % single solution
    muEq=subssolve(muEq{1},phi_variables,phi_solution,mpnb);
    
end

if numerical
  % convert to doubles  
  for i=length(muEq):-1:1
      try
          muEq{i}=double(muEq{i});
          if numerical<0    
              if any(any((muEq{i}<0) | (abs(imag(muEq{i}))>1e-6) ))
                  muEq(i)=[];      
                  if exist('phiEq','var')  
                      phiEq(i)=[];
                  end	
                  if exist('xiEq','var')  
                      xiEq(i)=[];
                  end	
              end      
          end      
          if exist('phiEq','var')  
              for i=1:length(phiEq)
                  phiEq{i}=double(phiEq{i});    
              end  
          end    
      catch
          warning('evaluation of numerical solution %d does not return double',i);
      end
  end  
end  

fprintf('finished %.2fsec\n',etime(clock,t0));

%delete(mpnb)  % should be used, but apparently creates errors
