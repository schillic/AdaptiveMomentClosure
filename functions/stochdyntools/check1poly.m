function rc=check1poly(symSpecies,symParameters,expr,str,mpnb)

rc=1;

try
  coef=monomialindexes(sym(char(expr)),symSpecies,mpnb);
catch
  fprintf('checkPolynomial: the following %s is not a polynomial on the species\n\t %s: ''%s''\n',str,str,char(expr));
  fprintf('\t species: ''%s''\n',findsym(symSpecies))
  rc=0;
  return 
end

coefsym=sym(['[',findsym(coef),']']);
for i=1:length(coefsym)
  if all(symParameters~=coefsym(i))
    fprintf('checkPolynomial: unknown coefficient ''%s'' in the polynomial ''%s''\n',char(coefsym(i)),char(expr));
    rc=0;
  end
end