function out=subssolve(in,variables,solution,mpnb)
% out=subsolve(in,variables,solution)
%
% Takes the output of the Matlab command solve and uses it to
% eliminate the variables solved in the expression given
%
% Inputs
% ------
%
% in: matrix of symbolic expressions where replacements are needed
%
% variables: cell array of variables to be replaced
%            (as provided to the solve() command) 
%
% solution: output of the solve() command (single output)
%
% Output
% ------
%
% out: cell array of matrices of symbolic expressions with
%      replacements. One cell per solution.

if nargin<4
    mpnb=symengine;
end

%variables
%size(variables)

if length(variables)==1
  for i=1:length(solution)
    out{i}=mysubs(in,variables,solution(i),0,mpnb);
  end
  return  
end  

%in

for j=1:length(variables)
  thisvariable=getfield(solution,char(variables{j}));
  for i=1:length(thisvariable)
    if j==1    
        out{i}=mysubs(in,variables(j),thisvariable(i),0,mpnb);
    else
      out{i}=mysubs(out{i},variables(j),thisvariable(i),0,mpnb);
    end      
  end
end  

%out{1}

if nargin<4
    %delete(mpnb)  % should be used, but apparently creates errors
end
