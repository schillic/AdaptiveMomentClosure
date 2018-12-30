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

