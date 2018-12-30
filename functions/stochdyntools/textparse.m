function [rawLines,statements,statementLines]=textparse(filename);

fid=fopen(filename);

if fid<0
  error('error opening file ''%s''\n',filename)  
end

rawLines={};
statements={};
statementLines=[];
while 1
  tline = fgetl(fid);
  
  if ~ischar(tline), break, end
  
  rawLines(end+1)=cellstr(tline);
  
  %% remove comments
  k=find(tline=='%');
  if ~isempty(k)
    tline=tline(1:min(k));
  end
    
  %% break into statements (removing leading and traling white spaces)
  k=regexp(tline,'\s*+([^;:]*[^\s])\s*[;:]','tokens');
  for i=1:length(k)
    statements(end+1)=k{i};
    statementLines(end+1,1)=length(rawLines);
  end
end

fclose(fid);

