% Add injection losses to the squeezed field
% lambda_in is defined as ifo.Squeezer.InjectionLoss

function Mout = sqzInjectionLoss(Min, L)
  
  eye2 = eye(size(Min,1), size(Min,2));
  Meye = eye2(:, :, ones(size(Min, 3), 1));

  Mout = [Min .* sqrt(1 - L), Meye .* sqrt(L)];

end
