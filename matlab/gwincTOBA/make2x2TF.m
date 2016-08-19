% Create a transfer matrix with 2x2xnF
%   The vectors must all have nF elements
function M3 = make2x2TF(A11, A12, A21, A22)

  nF = max([numel(A11), numel(A12), numel(A21),numel(A22)]);

  % if any input is just a number, expand it
  if numel(A11) == 1
    A11 = A11 .* ones(1, nF);
  end
  if numel(A12) == 1
    A12 = A12 .* ones(1, nF);
  end
  if numel(A21) == 1
    A21 = A21 .* ones(1, nF);
  end
  if numel(A22) == 1
    A22 = A22 .* ones(1, nF);
  end
  
  % check to be sure that they all match now
  if any([numel(A11), numel(A12), numel(A21),numel(A22)] ~= nF)
     error('Input vector length mismatch.')
  end
  
  % build output matrix
  M3 = zeros(2,2,nF);
  for k = 1 : nF
    M3(:,:,k) = [A11(k) A12(k); A21(k) A22(k)];
  end
end
