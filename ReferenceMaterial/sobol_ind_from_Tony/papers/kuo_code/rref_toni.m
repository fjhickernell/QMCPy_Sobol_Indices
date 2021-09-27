function rrA = rref_toni(A, zeroThresh)

%  On entry: A is a square matrix,
%            b is a column vector of the same dimension,
%            display is 1 if step-by-step display desired, 0 otherwise.
%  On exit:  U and f are the result of applying Gaussian elimination
%            with partial pivoting to the system Ax=b.

[m, n] = size(A);
if (nargin == 1)
  zeroThresh = n*eps;
end

currRow = 1;
for j=1:n
  if (currRow <= m)
	[cmax,relPos] = max(abs(A(currRow:m,j)));  % index of max element in col j.
	temp = A(currRow,:);                     % Exchange rows j and j+jmax-1.
	A(currRow,:) = A(currRow + relPos - 1,:);
	A(currRow + relPos - 1,:) = temp;

	if (abs(A(currRow,j)) > zeroThresh)
	  A(currRow,:) = A(currRow,:)/A(currRow,j);
	  A(currRow,1:j-1) = zeros(1,j-1);
	  for i=1:currRow - 1                 % Eliminate in column j using currRow.
		mult = A(i,j);
		A(i,j:n) = mod(A(i,j:n) - mult*A(currRow,j:n),2);
	  end
	  for i=currRow + 1:m                 % Eliminate in column j using currRow.
		mult = A(i,j);
		A(i,j:n) = mod(A(i,j:n) - mult*A(currRow,j:n),2);
	  end
	  A([1:currRow-1 currRow+1:m],j) = zeros(m-1,1);
	  currRow = currRow + 1;
	else
	  A(currRow:m,j) = 0;
    end
  end
end
rrA = A;
end
