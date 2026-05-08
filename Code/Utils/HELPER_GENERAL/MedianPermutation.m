function [Median, MedPartition, Relations, Dom] = ...
    MedianPermutation(PermProfile, Weights, SingleMed)
  % Median permutation (order, ranking) from a permutation set (profile).
  %
  % function [Median, MedPartition, Relations, Dom] = ...
  %     MedianPermutation(PermProfile, Weights, SingleMed)
  %
  % Inputs:
  %   PermProfile: Matrix containing the set of N orders (rows) of M
  %     ranked items (columns) each for which to find one or all median
  %     orders.  The first item in an order is the most preferred.
  %   Weights [default = ones]: Vector of length N, to give different
  %     weights to each order.  The weights are treated as multiplicities
  %     of corresponding orders.  Thus they are multiplied in counts of
  %     relations and sum of distances.
  %   SingleMed [default = false]: Whether to return a single median
  %     permutation, even when more than one exist.
  %
  % Outputs: 
  %   Median: Matrix where each row is a median order of the provided
  %     profile.
  %   MedPartition: When more than one median exists, only a subset of
  %     items may exchange position.  This cell array gives a single
  %     "median-ordered partition" that represents all medians.  E.g. if
  %     there are three medians: [1,2,3,4,5], [1,3,4,2,5], [1,4,2,3,5],
  %     then MedPartition = {[1], [2, 3, 4], [5]}.
  %   Relations: Square matrix of size N, where element (e,f) indicates
  %     whether item e is preferred more often (+1), equally (0), or less
  %     often (-1) than item f among all orders in the given profile.
  %   Dom: Square symmetrical logical matrix of "dominating" relations,
  %     where an item is always preferred over another (or the opposite) in
  %     the given profile.
  %
  % Theory: A set of preference orders can be aggregated to obtain a single
  % representative order, sometimes called a consensus order.  One such
  % order is the median order, which minimizes the sum of distances between
  % itself and each order in the set.  This order, also called Kemeny, or
  % Kemeny-Young order, satisfies desirable properties such as the extended
  % Condorcet criterion [Truchon 1998]: if two subsets of scenarios are
  % such that when comparing one scenario from each subset a strict
  % majority of participants always prefers the scenario from the same
  % subset, then the consensus order should respect that preference, i.e.
  % have each scenario from the "majority preferred subset" ranked before
  % the ones from the other subset.  It can thus be considered a best
  % compromise order, taking each pairwise scenario preferences into
  % account.  It can also be derived as the maximum likelihood order if the
  % data is considered as a noisy sampling of a "correct" order.  Note that
  % a set of orders can have more than one median, each having the same sum
  % of distances.
  %
  % References: [Ali 2012, Experiments with Kemeny ranking: What works
  % when?], [Betzler 2009, How Similarity Helps to Efficiently Compute
  % Kemeny Rankings], [Bredereck 2009, Fixed-Parameter Algorithms for
  % Computing Kemeny Scores - Theory and Practice], [Truchon 1998, An
  % Extension of the Condorcet Criterion and Kemeny Orders]
  %
  % Marc Lalancette, 2015-10-04
  
  % Could there always be a descending or equal only path from any
  % permutation to the median? That would mean a better branching algorithm
  % than trying every permutation.  Seems likely though maybe not since
  % it's an NP-hard problem?
  
  %   % Increasing this parameter may improve performance by not partitioning
  %   % very small subsets.  Keep at 1 minimum. [NOT IMPLEMENTED]
  %   SubsetTargetSize = 1;
  %   if SubsetTargetSize < 1
  %     SubsetTargetSize = 1;
  %   end
  
  % Option to return only one median permutation, even when more than one
  % exist.
  if nargin < 3 || isempty(SingleMed) || ~SingleMed
    SingleMed = false;
  else
    SingleMed = true;
  end
  
  % nP permutations (rankings, orders, preference lists, votes), nE elements (ranked items, candidates)
  [nP, nE] = size(PermProfile); 
  [~, Positions] = sort(PermProfile, 2);
  
  % Treat weights as multiplicities of corresponding orders.  Thus they are
  % multiplied in counts of relations and sum of distances.
  if ~exist('Weights', 'var') || isempty(Weights)
    Weights = ones(nP, 1);
  else
    % Ensure they sum up to nP.
    Weights = nP / sum(Weights) * Weights;
    % Ensure they are a column.
    if ~iscolumn(Weights)
      Weights = Weights';
    end
  end
  
  % Pairwise comparison of positions.
  PP = zeros(nE);
  for e = 1:nE-1
    for f = e+1:nE
      % Count of number of times pos(e) > pos(f) - count of pos(e) < pos(f).
      % == 0 is possible with even number of permutations.
      PP(e, f) = nP - 2 * sum(Weights .* (Positions(:, e) < Positions(:, f))); 
    end
  end
  PP = PP - PP';
  
  % Element (row)'s position is greater than element (column)'s.  It comes
  % after.  This is used in subfunction that calculates distances of all
  % possible candidates.  More efficient to do it once globally.
  
  % The "precedence matrix" is sufficient for determining the Kemeny order,
  % so there's likely a way to not have to loop over the profile. [e.g. Ali
  % 2012]  Could After be summed here already over nP (with weights)?

  % Weights are multiplied later.
  After = false(nE, nE, nP);
  for e = 1:nE-1
    for f = e+1:nE
      After(e, f, :) = Positions(:, e) > Positions(:, f);
      After(f, e, :) = ~After(e, f, :);
    end
  end
  %   After = double(After);
  % probably faster to sum logical than double?
  
  % Identify non-dirty (dominating) relations, those that are always the
  % same in each permutation.  They are maintained in the median. [Betzler
  % 2009, apparently follows from the extended Condorcet criterion.]
  Dom = abs(PP) == nP;
  Relations = sign(PP); % sign(0) = 0.
  clear PP % PermProfile
  
  % Order will first be partitioned maximally into subsets that satisfy the
  % extended Condorcet criterion [Truchon 1998]: they are in a consistent
  % (no cycles) strict majority order, which is maintained in the median.
  % Working order.
  Median = 1:nE;
  % Indices of start of subsets (end before next start, or last element).
  Partition = [1, nE+1];
  % Subsets that are already maximally partitioned.  Define these as
  % non-dirty subsets: each member only has non-dirty relations with all
  % other elements not in the subset.  (Different from [Bredereck 2009
  % thesis], where he defines a dirty set as having a connected graph of
  % dirty relations, so looking at relations within the set, whereas I
  % consider external relations.)
  Clean = false;
  
  % Find non-dirty elements (not part of any dirty relations thus always
  % same position with same subsets on each side) position and split.
  % These are the singletons in the desired partition.
  CleanE = find(sum(Dom, 1) == nE-1); % linear indices
  for e = CleanE
    % Get subset.
    [iP, L] = GetSubset(e);
    % Find subsets on each side.
    [Smaller, Equal, Greater] = RelationSets(L, e); % Equal == e.
    % Split
    SplitG = ~isempty(Greater);
    SplitS = ~isempty(Smaller);
    if SplitG || SplitS
      Median(Partition(iP):Partition(iP+1)-1) = [Smaller, Equal, Greater];
      if SplitG && SplitS
        Partition = [Partition(1:iP), Partition(iP) + ...
          numel(Smaller) + [0, 1], Partition(iP+1:end)];
        Clean = [Clean(1:iP-1), [false, true, false], Clean(iP+1:end)];
      elseif SplitG
        Partition = [Partition(1:iP), Partition(iP) + 1, ...
          Partition(iP+1:end)];
        Clean = [Clean(1:iP-1), [true, false], Clean(iP+1:end)];
      else % if SplitS
        Partition = [Partition(1:iP), Partition(iP+1) - 1, ... % or iP) + numel(Smaller)
          Partition(iP+1:end)];
        Clean = [Clean(1:iP-1), [false, true], Clean(iP+1:end)];
      end
    else % if ~Clean(iP) % Always, otherwise would mean did same one twice.
      % Singleton already.
      Clean(iP) = true;
    end
  end
  
  % For each subset that is not known to be "clean" (finest partition
  % possible), find the clean subset for the first element and split
  % accordingly.
  while any(~Clean)
    iP = find(~Clean, 1, 'first');
    % Mark clean if only 1 element.
    if Partition(iP+1) - Partition(iP) == 1
      Clean(iP) = true;
      continue;
    end
    
    L = Median(Partition(iP):Partition(iP+1)-1);
    
    % Finest subset for element c, that will be in the final partition.
    [SmallerEq, ~, GreaterEq] = FullRelationSets(L, L(1));
    % The Equal subset of FullRelationSets is the connected graph of
    % equalities only.  Below, Equal will contain all cycles involving
    % those elements.  Greater and Smaller are in strict majority relations
    % with all elements of the other two sets.
    Equal = intersect(GreaterEq, SmallerEq);
    Greater = setdiff(GreaterEq, Equal);
    Smaller = setdiff(SmallerEq, Equal);

    SplitG = ~isempty(Greater);
    SplitS = ~isempty(Smaller);
    if SplitG || SplitS
      Median(Partition(iP):Partition(iP+1)-1) = [Smaller, Equal, Greater];
      if SplitG && SplitS
        Partition = [Partition(1:iP), Partition(iP) + ...
          numel(Smaller) + [0, numel(Equal)], Partition(iP+1:end)];
        Clean = [Clean(1:iP-1), [false, true, false], Clean(iP+1:end)];
      elseif SplitG
        Partition = [Partition(1:iP), Partition(iP) + numel(Equal), ...
          Partition(iP+1:end)]; % + numel(Smaller) == 0
        Clean = [Clean(1:iP-1), [true, false], Clean(iP+1:end)];
      else % if SplitS
        Partition = [Partition(1:iP), Partition(iP) + numel(Smaller), ...
          Partition(iP+1:end)];
        Clean = [Clean(1:iP-1), [false, true], Clean(iP+1:end)];
      end
    else % if ~Clean(iP) % Always, wouldn't have tried it otherwise.
      Clean(iP) = true;
    end
  end
  
  % ---------------------------------------------------------------------
  % Partitioning complete.  Process each subset separately.
  nPart = numel(Partition) - 1;
  % Each subset can have multiple median (partial) orders, so prepare cells
  % to store each unknown-sized list of orders.
  MedPartition = cell(1, nPart);
  MedPartCount = zeros(1, nPart);
  PartSize = Partition(2:end) - Partition(1:end-1);
  MaxPartSize = max(PartSize);
  if MaxPartSize > 10
    error('Max partition size is %d, probably too large for this function.', ...
      MaxPartSize);
  end
  for iP = 1:nPart
    if PartSize(iP) == 1
      MedPartition{iP} = Median(Partition(iP));
      MedPartCount(iP) = 1;
      continue;
    else
      % Need to test all possible permutations (that respect dominating
      % relations).
      
      % Actually, here the permutations also need to respect consecutive
      % majority relations.  There are many more of these, so it seems it
      % would be much better to build them directly with a recursive
      % function instead of using perms and pruning.  But no time right
      % now.
      
      L = Median(Partition(iP):Partition(iP+1)-1);
      % Get all permutations.
      SubPerm = perms(1:PartSize(iP));

      if any(any(Relations(L, L)))
        % Check if there are dominating relations within the subset to
        % filter unneeded permutations.
        DRel = Dom(L, L) .* Relations(L, L);
        if any(DRel(:))
          % Sorting the full set of permutations is expensive, but hadn't
          % thought of another way, unless we ignore dominating relations
          % (but now see comment above).
          [~, SubPos] = sort(SubPerm, 2);
          % Only need half the relations (since they are anti-symmetrical).
          [iLess, iGreat] = find(DRel < 0);
          Reject = false(size(SubPos, 1), 1);
          for i = 1:numel(iLess)
            Reject = Reject | SubPos(:, iLess(i)) > SubPos(:, iGreat(i));
          end
          clear SubPos
          SubPerm(Reject, :) = [];
        end
        
        % Apply order to subset, now we need the order, and the positions
        % of the profile permutations.
        SubPerm = L(SubPerm);
        % Loop on profile.  (Looping on SubPerm was incredibly longer.) 
        % (However see comment near After definition above.)
        D = PartialDistanceProfileLoop(SubPerm);
        % SubScore = min(D);
        SubPerm = SubPerm(D == min(D), :);
        clear D
      else % All equal.
        % All permutations valid.
        SubPerm = L(SubPerm);
      end
      % Store valid orders.
      MedPartition{iP} = SubPerm;
      MedPartCount(iP) = size(SubPerm, 1);
    end
  end % Partition loop
  
  % Organize output list
  MedPartMult = cumprod(MedPartCount)';
  nMed = MedPartMult(end);
  MedPartMult = [[1; MedPartMult(1:end-1)], nMed ./ MedPartMult];
  Median = zeros(nMed, nE);
  for iP = 1:nPart
    for r = 1:MedPartMult(iP, 2)
      for m = 1:MedPartCount(iP)
        Median(((r-1)*MedPartCount(iP) +(m-1))*MedPartMult(iP, 1) + (1:MedPartMult(iP, 1)), ...
          Partition(iP):Partition(iP+1)-1) = ...
          MedPartition{iP}(ones(MedPartMult(iP, 1), 1) * m, :);
      end
    end
  end
  
  if SingleMed
    Median = Median(1, :);
  end
  
  
  % ---------------------------- subfunctions that use global variables.
  function [iP, L] = GetSubset(c)
    iM = find(Median == c, 1, 'first');
    iP = find(Partition <= iM, 1, 'last');
    if nargout > 1
      L = Median(Partition(iP):Partition(iP+1)-1);
    end
  end
  
  function [Smaller, Equal, Greater] = RelationSets(L, c)
    Greater = L(Relations(c, L) < 0);
    Smaller = L(Relations(c, L) > 0);
    Equal = L(Relations(c, L) == 0);
  end
  
  % Same as above, but iteratively on subset members.
  function [SmallerEq, Equal, GreaterEq] = FullRelationSets(L, c)
    % Find '=' set first.  (These are always in the same partition.)
    Equal = c;
    Add = setdiff(L(any(Relations(Equal, L) == 0, 1)), Equal);
    while ~isempty(Add)
      Equal = [Equal, Add]; %#ok<*AGROW>
      Add = setdiff(L(any(Relations(Equal, L) == 0, 1)), Equal);
    end
    GreaterEq = Equal;
    Add = setdiff(L(any(Relations(GreaterEq, L) <= 0, 1)), GreaterEq);
    while ~isempty(Add)
      GreaterEq = [GreaterEq, Add];
      Add = setdiff(L(any(Relations(GreaterEq, L) <= 0, 1)), GreaterEq);
    end
    SmallerEq = Equal;
    Add = setdiff(L(any(Relations(SmallerEq, L) >= 0, 1)), SmallerEq);
    while ~isempty(Add)
      SmallerEq = [SmallerEq, Add];
      Add = setdiff(L(any(Relations(SmallerEq, L) >= 0, 1)), SmallerEq);
    end
  end
    
  % Same as PartialDistance, but loop over profile permutations and sum for
  % all reference permutations simultaneously.  Much faster, but can no
  % longer stop as soon as D >= MaxD.
  function D = PartialDistanceProfileLoop(RefPerm)
    % Count number of distinct element pairs in wrong order, equivalent to
    % counting the number of swaps (interchanging two adjacent elements)
    % required to (partially) sort the permutation.
    [nPSub, nESub] = size(RefPerm);
    D = zeros(nPSub, 1);
    for pp = 1:nP
      ppD = zeros(nPSub, 1);
      ppAfter = After(:, :, pp);
      
      for ee = 1:nESub-1
        % Use linear indices into position comparison matrix. 
        for ff = ee+1:nESub
          ppD = ppD + ppAfter(RefPerm(:, ee) + nE .* (RefPerm(:, ff) - 1));
        end
      end
      D = D + Weights(pp) * ppD;
    end
  end
    
end



















