function Distances = PermutationDistance(Perm, ReferencePerm)
  % Kendall tau (bubble-sort, swap) distance between permutations (orders, rankings).
  % 
  % Distances = PermutationDistance(Perm, ReferencePerm)
  %
  % Calculate the number of adjacent element swaps to go from each
  % permutation (row) in the Perm array to one or more permutation in the
  % ReferencePerm array. If more than one line in ReferencePerm, the
  % average distance to the references is given for each permutation in
  % Perm.
  %
  % Theory: A widely used distance metric on the space of permutations is
  % the Kendall tau distance, or bubble-sort or swap distance.  It
  % corresponds to the number of pairwise disagreements between the two
  % considered orders.  The distance also corresponds to the number of
  % times two adjacent scenarios must swap positions to make one order
  % match the other.
  %
  % Marc Lalancette, 2014-02-10

  % Could be various ways to loop optimally depending on which of nP, nR,
  % nE are larger.
  
  [nP, nE] = size(Perm); % nP permutations (participants), nE elements (ranked items)
  nR = size(ReferencePerm, 1);
    
  Distances = zeros(nP, 1);
  for iR = 1:nR
    % Apply inverse reference permutation to individual permutations.
    [~, InverseReference] = sort(ReferencePerm(iR, :), 2);
    InvRefPerm = InverseReference(Perm);
    
    % Calculate distance between individual permutations and identity:
    % sortedness (count number of distinct element pairs in wrong order),
    % equivalent to counting the number of swaps (interchanging two
    % adjacent elements) required to sort the permutation.
    for e = 1:nE-1
      for f = e+1:nE
        Distances = Distances + (InvRefPerm(:, e) > InvRefPerm(:, f));
      end
    end
  end
  Distances = Distances / nR;
  
end
