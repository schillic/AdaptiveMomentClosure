function latex2clip(expr,msg,stop)

str=latex(expr);
str=regexprep(str,'{\\it (mu|nu|rho)\\_([a-zA-Z0-9])}','\\$1_$2');
str=regexprep(str,'{\\it (mu|nu|rho)\\_([a-zA-Z0-9]+)}','\\$1_{$2}');
str=regexprep(str,'{\\it ([a-zA-Z0-9]+)\\_([a-zA-Z0-9])}','$1_$2','preservecase');
str=regexprep(str,'{\\it ([a-zA-Z0-9]+)\\_([a-zA-Z0-9]+)}','$1_{$2}','preservecase');
str=regexprep(str,'{\\it ([a-zA-Z0-9]+)}','$1','preservecase');
str=regexprep(str,'{([0-9])}','$1');

clipboard('copy',str);
if stop
  fprintf('%s copied to clipboard <paused>\n',char(msg))
  pause  
end  