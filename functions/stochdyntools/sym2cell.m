function cel=sym2cell(arr)
% cel=sym2cell(arr)
% 
% Converts a symbolic array to an array of cells of the same size
% 
% For numerical arrays this can be done simply by cel={arr},
% but this option does not seem to work for symbolic arrays
    
cel=mat2cell(arr,ones(size(arr,1),1),ones(size(arr,2),1));
