function [nK, ProbMemb, Med, Spread, ProbClust, LogLikes, AICcs] = ...
    MallowsClustering(Perm, nRep, nKs, ...
    LogLikeTol, IterTol, Display, Verbose)
  % Parametric expectation maximization clustering for permutations using Mallows' phi distribution.
  %
  % [nK, ProbMemb, Med, Spread, ProbClust, LogLikes, AICcs] = ...
  %     MallowsClustering(Perm, nRep, nKs, LogLikeTol, IterTol, Display, Verbose)
  %
  % Inputs: 
  %   Perm: Matrix containing the set of N orders (rows) of M ranked items
  %     (columns) each.  The first item in an order is the most preferred.
  %   nRep: Number of times to repeat with different initial conditions,
  %     i.e. cluster centers.
  %   nKs: Number of clusters to try, as a vector, e.g. 1:5.  Each will be
  %     repeated nRep times and the best of each compared to find the
  %     overall best cluster solution.
  %   LogLikeTol [default = 1e-10]: Stopping criterion on the change of
  %     log-likelihood.
  %   IterTol [default = 200]: Max number of iterations of the EM algorithm
  %     (for each repetition).
  %   Display [default = 0]: Option to plot the clustering process and/or
  %     result, where each order is a circle colored according to it's
  %     "preferred" cluster.  The circles are positioned according to the
  %     first two PCA axes of the profile.  Accepted values: 0 no display,
  %     1 final result only, 2 each repetition, 3 each EM iteration.
  %   Verbose [default = 0]: Print some information on progress in command
  %     window.
  %
  % Outputs: 
  %   nK: Number of clusters of the best found solution.  Even if some
  %     cluster numbers are skipped (e.g. nKs = [1, 3, 5]), the last two
  %     outputs will be indexed by nK, thus have a size of max(nKs).
  %   ProbMemb: Cluster membership probabilities of each order in each
  %     cluster ("soft" clustering).
  %   Med: Median (central) order of each cluster.
  %   Spread: Dispersion parameter of the Mallows' phi distribution for 
  %     each cluster.
  %   ProbClust: Probability of a random order (according to the best 
  %     model) being a member of each cluster.
  %   LogLikes: Log-likelihood of the model at the end of each repetition
  %     of the EM algorithm, for each number of cluster (size= [max(nKs,
  %     nRep]).
  %   AICcs: Akaike information criterion of each repetition (same size as
  %     LogLikes).  This is used to compare the best models of each number
  %     of cluster, to select the best overall model.  
  %
  % Theory and references: Mallows' phi distribution is an exponential
  % probability distribution in the space of permutations, based on the
  % Kendall tau distance, and having as parameters a central order, and a
  % dispersion parameter.  We use the median order as the central
  % parameter.  A single distribution may not be appropriate to represent
  % the data, thus this function compares mixture models with various
  % numbers of components or clusters. The best model for a given number of
  % clusters is found using an expectation maximization (EM) procedure, as
  % in [Murphy 2003, Mixtures of distance-based models for ranking data] or
  % [Busse 2007, Cluster Analysis of Heterogeneous Rank Data].  This
  % process starts with an initial choice for the parameters and then
  % iteratively computes the cluster membership probabilities of each order
  % (expectation step), and the parameters that maximize the expectation
  % formula (maximization step).
  % AIC is a bias-corrected estimator of expected relative information loss
  % when replacing the data with the model.  For a discussion of AIC in
  % model selection, in particular comparing it to BIC, see [Burnham 2004,
  % Multimodel Inference, Understanding AIC and BIC in Model Selection].
  %
  % Marc Lalancette, 2015-10-11
  
  % Saw at least in 2 places a stopping criterion of diff(log-likelihood) <
  % 1e-10, but seems a relative tolerance would make more sense.
  
  % [I was hoping to have a slight difference with Mallows' model in that
  % we'd use the average distance to all median permutations instead of a
  % single modal permutation. In practice, however, since it is a weighted
  % median, I think it is very unlikely that there would be multiple
  % medians in most cases.  And it makes the partition function not
  % simplifiable (at least not obviously) so we'd have to calculate it over
  % the symmetrical group...]
  
  if nargin < 7 || isempty(Verbose)
    Verbose = 0;
  end
  if nargin < 6 || isempty(Display)
    Display = 0;
  end
  if nargin < 5 || isempty(IterTol)
    IterTol = 200;
  end
  if nargin < 4 || isempty(LogLikeTol)
    LogLikeTol = 1e-10;
  end
  
  nKmax = max(nKs(:));
  
  [nP, nE] = size(Perm);
  if nE < 2
    Display = 0;
  end
  if Display
    % Do PCA to plot 2 principal components.
    
    [~, Ranks] = sort(Perm, 2);
    g = 1/(nP - 1);
    RanksC = bsxfun(@minus, Ranks, mean(Ranks, 1));
    Cov = g * (RanksC' * RanksC);
    YTotalVar = trace(Cov);
    
    [~, PCA_LSqrt, PCA_V] = svd(sqrt(g) * RanksC, 'econ'); %
    PCA_L = (diag(PCA_LSqrt).^2)';
    PCA_VarProportion2 = sum(PCA_L(1:2)) / YTotalVar;
    PCA_Comp = RanksC * PCA_V; % Verified
    
    plot(PCA_Comp(:, 1), PCA_Comp(:, 2), 'ok');
    hold on
    % set(gca, 'XLim', [-8, 8], 'YLim', [-7, 7]);
    xlabel('1^{st} PCA component of scenario rankings');
    ylabel('2^{nd} PCA component');
    title( {'Scatter of permutations, as projected on', ...
      sprintf( 'first 2 PCA components (%d%% of variance).', ...
      round(100*PCA_VarProportion2) )} );
    
    Colors = colormap(lines(nKmax));
    H = zeros(nKmax, 1);
    for k = 1:nKmax
      H(k) = plot(-inf, -inf, 'o', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', Colors(k, :));
    end
    drawnow;
  end
  
  % Probabilies.
  %   ProbPerm = zeros(nP, 1); % P(X)
  ProbPermInClust = zeros(nP, nKmax); % P(X|Z)
  ProbMembership = zeros(nP, nKmax); % P(Z|X)
  % Model parameters.
  ProbClust = zeros(1, nKmax); % P(Z)
  Med = zeros(nKmax, nE);
  Spread = zeros(1, nKmax);
  % Returned results.
  BestProbMemb = zeros(nP, nKmax); % P(Z|X)
  BestProbClust = zeros(1, nKmax); % P(Z)
  BestMed = zeros(nKmax, nE);
  BestSpread = zeros(1, nKmax);
  %   BestProbPermInClust = zeros(nP, nKmax);
  BestAICc = inf;
  
  if nargout >= 6
    DoLogLikes = true;
    LogLikes = zeros(nKmax, nRep);
    AICcs = zeros(nKmax, nRep);
  else
    DoLogLikes = false;
  end
  
  % Loop over given list of number of clusters/components.
  for nK = nKs
    
    if Verbose > 1
      fprintf('\nk: %d\n  rep:    0', nK);
    end
    % Repeat with many different starting conditions.  The number of local
    % maxima found across repetitions increases greatly with more clusters.
    % But that could simply be due to our very limited initialization
    % choices with few clusters.  Perhaps should add a "noise" element to
    % the starting medians (say a few swaps).  Would also benefit from
    % adjusting nRep based on that, to say 100 times the number found so
    % far, to be confident of finding the global maximum.
    for r = 1:nRep
      if Verbose > 1
        fprintf('\b\b\b\b%4d', r);
      end
      if r > 1 && nK == 1
        % No need for repetitions with one cluster.
        break;
      end
      LogLikeDiff = LogLikeTol + abs(LogLikeTol); % Ensure we start above stopping criterion.
      PrevLogLike = -inf;
      Iter = 0;
      
      % Initialize by choosing k elements as cluster centers and a
      % spread of the inverse average distance to all elements.
      iP = randperm(nP, nK); % nK distinct integers in 1:nP.
      Med(1:nK, :) = Perm(iP, :);
      for k = 1:nK
        Spread(k) = 1 / mean( PermutationDistance(Perm, Perm(k, :)) );
      end
      % Uniform cluster probabilities.
      ProbClust(1:nK) = 1/nK;
      
      while Iter < IterTol && LogLikeDiff > LogLikeTol
        % (LogLikeDiff < 0 || ) % Not stopping if getting worse, mostly for debugging.  There is a warning if this happens, below.
        
        % -----------------------------------------------------
        %   E step
        % Calculate membership probabilities.
        
        % P(X|Z) Conditional probability of a data point given membership
        % in a cluster.
        for k = 1:nK
          % Single cluster pdf according to parametric model.
          e = exp(-Spread(k));
          ProbPermInClust(:, k) = ...
            e .^ PermutationDistance(Perm, Med(k, :)) ./ ...
            prod((1 - e .^ (1:nE)) ./ (1 - e));
          % With one random initialization, 1 cluster, large spead (on the
          % order of average distance to cluster members) means the
          % numerator for a random permutation is on the order of 0.01 on
          % average.  15! = 1.3e12, and the denominator was verified to
          % give 1.4e10.  All consistent.  Numerator max is 1 (at the
          % median), average is 0.3.  So extremely small probability of
          % obtaining any of Perm with "large" spread is correct.
          %             ClustProbDist(Perm, Med(k, :), Spread(k), nE);
        end
        
        % P(X)
        ProbPerm = sum( ...
          bsxfun(@times, ProbPermInClust(:, 1:nK), ProbClust(1:nK)), 2 );
        
        % Likelihood of the model (parameters) given the data = probability
        % of obtaining the data given the model P(X) = probability with
        % cluster possibilities marginalized out = sum cluster
        % possibilities P(X,Z) = sum P(Z)P(X|Z)
        %  correct if cluster possibilities partition all of X, i.e. no
        %  possibility of x not being in a cluster.
        
        % This is the full log likelihood with the current estimate.  This
        % is what must increase.
        LogLike = sum(log(ProbPerm), 1);
        LogLikeDiff = LogLike - PrevLogLike;
        % Check that algorithm actually improves at each iteration.
        %         if LogLikeDiff < 0
        %           warning('Not improving; nK %d, iter %d.', nK, Iter);
        %         end
        
        %  P(Zi|Xj) = P(Xj|Zi)P(Zi)/P(Xj); P(Xj) = sum_i(P(Xj|Zi)P(Zi))
        % Can think of division by P(Xj), part of Bayes rule, as
        % normalizing such that each permutation has probability 1 of being
        % a member of a cluster.
        ProbMembership(:, 1:nK) = bsxfun(@rdivide, ...
          bsxfun(@times, ProbPermInClust(:, 1:nK), ProbClust(1:nK)), ...
          ProbPerm);

        % Check division by zero.
        %         if any(isnan(ProbMembership(:, 1:nK)))
        %           error('NaN detected, division by zero should not happen if properly initialized.');
        %         end
        
        % Update plot.
        if Display > 2
          Membership(:, 1:nK) = bsxfun(@eq, ...
            ProbMembership(:, 1:nK), max(ProbMembership(:, 1:nK), [], 2));
          for k = 1:nK
            set(H(k), 'xdata', PCA_Comp(Membership(:, k), 1), ...
              'ydata', PCA_Comp(Membership(:, k), 2));
          end
          drawnow;
          %         pause(0.5);
        end
        
        % -----------------------------------------------------
        %   M step
        % Calculate parameters that maximize or at least increase (for GEM)
        % the expectation formula.  Our parameters are the mixing factors =
        % cluster probabilities p(Z) and the median and spread for the
        % cluster distance probability density functions from which we get
        % p(X|Z).
        
        % Calculate mixing parameters (ProbK)
        %  P(Zi) = 1/n sum_j(P(Zi|Xj))
        % Probability of a random element being a member of each cluster.
        ProbClust(1:nK) = sum(ProbMembership(:, 1:nK), 1) ./ nP;

        for k = 1:nK
          % New spread.  Need to compute numerically.  Using distances from
          % previous iteration median.
          X = ProbMembership(:, k)' * PermutationDistance(Perm, Med(k, :)) ...
            ./ (nP * ProbClust(k));
          ClustSpreadFun = @(S) sum( 1/(exp(S) - 1) - (1:nE)./(exp((1:nE) .* S) - 1) ) - X;
          Spread(k) = fzero(ClustSpreadFun, Spread(k));

          % New median.
          Med(k, :) = MedianPermutation(Perm, ProbMembership(:, k), true);
        end
        
        % -----------------------------------------------------
        % Set for next iteration.
        Iter = Iter + 1;
        PrevLogLike = LogLike;

      end % Iterations for one repetition.
      
      if Display > 1
        Membership(:, 1:nK) = bsxfun(@eq, ...
          ProbMembership(:, 1:nK), max(ProbMembership(:, 1:nK), [], 2));
        for k = 1:nK
          set(H(k), 'xdata', PCA_Comp(Membership(:, k), 1), ...
            'ydata', PCA_Comp(Membership(:, k), 2));
        end
      end
      drawnow; % May make Matlab more responsive while running this?
      %         pause; % (0.5)
      
      if Iter >= IterTol
        warning('Clustering did not converge within %d iterations. Last log likelihood diff = %g', ...
          IterTol, LogLikeDiff);
      elseif Verbose == 1
        fprintf('Clustering converged in %d iterations. \nLast log likelihood diff = %g\n', ...
          Iter, LogLikeDiff);
      end
      
      %  -----------------------------------------------------
      % Compare to previous best and keep best.  Use Akaike criterion.
      % [Burnham 2004] for discussion of AIC vs BIC.  [Murphy 2003]
      % compares BIC and ICL.
      
      % Number of freely estimatable model parameters: median and spread
      % for each cluster, plus nK-1 for relative cluster probabilities (-1
      % because constrained to sum to 1).
      nFree = 3 * nK - 1; 
      AICc = -2 * LogLike + 2 * nFree * (1 + (nFree + 1)/(nP - nFree - 1));
      if DoLogLikes
        LogLikes(nK, r) = LogLike;
        AICcs(nK, r) = AICc;
      end
      
      if AICc < BestAICc
        BestAICc = AICc;
        BestnK = nK;
        BestProbMemb(:, 1:nK) = ProbMembership(:, 1:nK);
        BestProbClust(1:nK) = ProbClust(1:nK);
        BestMed(1:nK, :) = Med(1:nK, :);
        BestSpread(1:nK) = Spread(1:nK);
        %         BestProbPermInClust(:) = 0;
        %         BestProbPermInClust(:, 1:nK) = ProbPermInClust(:, 1:nK);
      end
      
    end % Repetitions
  end % nK loop
  
  % Display best case.
  if Display
    Membership(:, 1:nK) = bsxfun(@eq, ...
      BestProbMemb(:, 1:nK), max(BestProbMemb(:, 1:nK), [], 2));
    for k = 1:nK
      set(H(k), 'xdata', PCA_Comp(Membership(:, k), 1), ...
        'ydata', PCA_Comp(Membership(:, k), 2));
    end
    drawnow;
  end
  
  % Return best case parameters.
  nK = BestnK;
  ProbMemb = BestProbMemb(:, 1:nK);
  Med = BestMed(1:nK, :);
  Spread = BestSpread(1:nK);
  ProbClust = BestProbClust(1:nK);
  
  if Verbose > 1
    fprintf('\nDone!\n');
  end
  
end

% ---------------------------------------------------------------------
% % Single cluster pdf according to parametric model.
% function p = ClustProbDist(Perm, Med, Spread, nE)
%   e = exp(-Spread);
%   p = e .^ PermutationDistance(Perm, Med) ./ ...
%     prod((1 - e .^ (1:nE)) ./ (1 - e));
% end
% 






