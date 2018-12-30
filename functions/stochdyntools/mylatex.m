function str=mylatex(expr,texrules)
% str=mylatex(expr,texrules)
%
% This function returns the LaTeX representation of the symbolic
% expression expr. 
% 
% Joao Hespanha's version of Matlab's latex() command, with several
% "beautifications." This function also allows the expression to be
% transformed by a set of replacement rules.
%
% Inputs
% ------
%
% expr : Vector of symbolic expressions
%
% texrules : cell array of strings with one row per
%            transformation rule and two columns.
%            . the first column contains a regular expression, following
%              the syntax of Matlab's regexp() 
%            . the second column contains the string to replace the
%              regular expression 
%
%            The transformation rules are applied sequentially and only
%            once. Replacements are not applied to strings already replaced.
%
% Outputs
% -------
%
% str : String with latex output 
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
% Created May 29, 2007 
%
% Modified May 20, 2009
% Removed \mathrm from variable names

if nargin==0
  str=sprintf('\\usepackage{amsmath,amsfonts}\n\\newcommand{\\matt}[1]{\\begin{bmatrix}#1\\end{bmatrix}}\n\\DeclareMathOperator{\\E}{E}');
  return;  
end

str=latex(expr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove \mathrm{ ... } for variable names 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str=regexprep(str,'\\mathrm\{(\w+)\}','$1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% use \matt{ ... } for matrices instead of {array}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tokens=regexp(str,'(.*)\s*\\left[(\[]\s*\\begin\s*\{array\}\{c*\}(.*)\\end\s*\{array\}\s*\\right[)\]]\s*(.*)','tokens','once');

if length(tokens)>3
  error('mylatex: expression appears to have matrices within matrices\n''%s''\n',str);
elseif  length(tokens)==3
   % remove   \noalign{\medskip} after \\
   [starts,ends]=regexp(tokens{2},'\\\\\\noalign\{\\medskip\}','start','end');
   for i=length(ends):-1:1
     tokens{2}=sprintf('%s\\\\\n%s',tokens{2}(1:starts(i)-1),tokens{2}(ends(i)+1:length(tokens{2})));
   end     
   
   str=sprintf('%s \\matt{ %s } %s',tokens{1},tokens{2},tokens{3});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove \,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str=regexprep(str,'\\,','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove {\it ... } for multi-letter variable names 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[starts,ends,tokens]=regexp(str,'\{\\it\s*([^}]+)\s*\}','start','end','tokens');
for i=length(ends):-1:1
  str=[str(1:starts(i)-1),' ',char(tokens{i}),' ',str(ends(i)+1:length(str))];
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove { ... } around variable names -- problematic because of \frac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[starts,ends,tokens]=regexp(str,'\{\s*(\w+)\s*\}','start','end','tokens');
%for i=length(ends):-1:1
%  str=[str(1:starts(i)-1),' ',char(tokens{i}),' ',str(ends(i)+1:length(str))];
%end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove \left( ... \right) around variable names 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[starts,ends,tokens]=regexp(str,'\\left\(\s*(\w+)\s*\\right\)','start','end','tokens');
for i=length(ends):-1:1
  str=[str(1:starts(i)-1),' ',char(tokens{i}),' ',str(ends(i)+1:length(str))];
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove \left( ... \right) around expressions without ()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[starts,ends,tokens]=regexp(str,'\\left\(\s*([^()]+)\s*\\right\)','start','end','tokens');
for i=length(ends):-1:1
  str=[str(1:starts(i)-1),' (',char(tokens{i}),') ',str(ends(i)+1:length(str))];
end     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% apply rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>1

  pos=[];

  % find all items to replace & replace then by zeros (not be found again)
  for k=1:size(texrules,1)
    [starts,ends]=regexp(str,texrules{k,1},'start','end');
    for i=length(ends):-1:1
      str(starts(i):ends(i))=0;
    end     
    pos=[pos,[starts;ends;k*ones(1,size(starts,2))]];
  end

  [dummy,k]=sort(pos(1,:),2,'descend');
  pos=pos(:,k);

  % now do all the replacements
  for i=1:size(pos,2)
    str=[str(1:pos(1,i)-1),texrules{pos(3,i),2},str(pos(2,i)+1:end)];   
  end    
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove double spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[starts,ends]=regexp(str,'  ','start','end');
while ~isempty(ends)
  for i=length(ends):-1:1
    str=[str(1:starts(i)-1),' ',str(ends(i)+1:length(str))];
  end     
  [starts,ends]=regexp(str,'  ','start','end');
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove spaces and the beginning & end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while str(1)==' '
  str=str(2:end);  
end  
while str(end)==' '
  str=str(1:end-1);  
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove spaces around { } ( ) + - ^ &
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[starts,ends]=regexp(str,' [\[\]{}()+^&\-]','start','end');
for i=length(ends):-1:1
  str=[str(1:starts(i)-1),str(ends(i):length(str))];
end     

[starts,ends]=regexp(str,'[\[\]{}()+^&\-] ','start','end');
for i=length(ends):-1:1
  str=[str(1:starts(i)),str(ends(i)+1:length(str))];
end     

[starts,ends]=regexp(str,' \\\\','start','end');
for i=length(ends):-1:1
  str=[str(1:starts(i)-1),str(ends(i)-1:length(str))];
end     

