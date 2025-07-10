function [datesCDS, survProbs, intensities] = bootstrapCDS_NOaccrual(datesDF, discounts, datesCDS, spreadsCDS, recovery)
% bootstrapCDS_NOaccrual computes survival probabilities and hazard rates 
% (intensities) for a set of CDS contracts, neglecting the accrual on the
% premium leg.
%
% INPUT:
%   datesDF    : Vector of dates corresponding to discount factors.
%   discounts  : Vector of discount factors.
%   datesCDS   : Vector of CDS expiry/maturity dates.
%   spreadsCDS : Vector of CDS spreads.
%   recovery   : Recovery rate for the CDS.
%
% OUTPUT:
%   datesCDS   : The same CDS expiry dates as input.
%   survProbs  : Computed survival probabilities at each CDS date.
%   intensities: Hazard rates (intensities) corresponding to each CDS date.

%% PRELIMINARY CALCULATIONS -----------------------------------------------
% Define the settlement date as the first discount factor date.
settlement = datesDF(1);

% Convert discount factor dates to year fractions from settlement using
% day count convention 3 (e.g., Act/365 or similar).
dates_zeroRates = yearfrac(settlement, datesDF, 3);

% Compute zero rates from the discount factors.
zRates = zeroRates(datesDF, discounts);

% Convert CDS expiry dates into year-frac relative to the settlement date.
dates_CDS_zrates = yearfrac(settlement, datesCDS, 3);

% Interpolate zero rates at the CDS expiry dates and convert them into 
% discount factors. The /100 assumes the zero rates are expressed in %.
discounts_CDS = exp(-interp1(dates_zeroRates, zRates, dates_CDS_zrates)...
                    /100 .* dates_CDS_zrates);

% Compute the time increments (delta_t) for each CDS maturity.
% The 'diff' of the array [0; yearfrac(...)] gives the time length of each 
% period using day count convention 6 (e.g., 30/360).
delta_t = diff([0; yearfrac(settlement, datesCDS, 6)]);

%% INITIALIZATION ---------------------------------------------------------
% Determine the number of CDS maturities/spreads.
N = length(spreadsCDS);

% Preallocate arrays for survival probabilities and hazard intensities.
survProbs = zeros(N, 1);
intensities = zeros(N, 1);

%% THEORY OF THE BOOTSTRAP PROCEDURE --------------------------------------
% For all year i in [1, N],
% [E]
% S_{iY} * SUM_{j=1}^{i} delta(t_{j-1},t_j) B(t_0,t_j) P(t_0,t_j)
% = 
% (1 - RR) * SUM_{j=1}^{i} B(t_0,t_j) [P(t_0,t_{j-1}) - P(t_0,t_j)]
% i.e., PREMIUM LEG for maturity i Y = PROTECTION LEG for maturity i Y.
% NB:  [P(t_0,t_{j-1}) - P(t_0,t_j)] = default happens in [t_{j-1},t_j].
% NB2: for i=1, only one term in each sum, simplify B(t_0,t_1), and get:
%      P(t_0,t_1) = (1 - RR) / [S_{1Y} * delta(t_0,t_1) + 1 - RR]
% NB3: for i>1, rewrite [E] the following way:
%      S_{iY} * RPV01(t_0,t_i)
%      =
%      (1 - RR) * SUM_{j=1}^{i} B(t_0,t_j) [P(t_0,t_{j-1}) - P(t_0,t_j)]
%      i.e.:
%      S_{iY} * RPV01(t_0,t_{i-1}) + S_{iY} * delta(t_{i-1},t_i) B(t_0,t_i)
%      P(t_0,t_i)
%      =
%      (1 - RR) * SUM_{j=1}^{i-1} B(t_0,t_j) [P(t_0,t_{j-1}) - P(t_0,t_j)]
%      + (1 - RR) * B(t_0,t_i) [P(t_0,t_{i-1}) - P(t_0,t_i)]
%      Then, isolate P(t_0,t_i) and you get:
%      P(t_0,t_i) = NUM / DEN where:
%      NUM = (1 - RR) * B(t_0,t_i) P(t_0,t_{i-1}) + (1 - RR) *
%            SUM_{j=1}^{i-1} B(t_0,t_j) [P(t_0,t_{j-1}) - P(t_0,t_j)]
%            - S_{iY} * RPV01(t_0,t_{i-1})
%      DEN = B(t_0,t_i) [(1 - RR) +  S_{iY} * delta(t_{i-1},t_i)]
% NB4: to obtain the intensity from the survival probability:
%      P(t_0,t_i) = exp(-SUM_{j=1}^{i}lambda(t_j)delta(t_{j-1},t_j)
%                 = P(t_0,t_{i-1}) * exp(-lambda(t_i)delta(t_{i-1},t_i))
%      i.e. lambda(t_i) = - 1 / delta(t_{i-1},t_i) * ln(P(t_0,t_{i}) /
%                                                       P(t_0,t_{i-1}))
%      i.e. lambda(t_i) = 1 / delta(t_{i-1},t_i) * ln(P(t_0,t_{i-1}) /
%                                                     P(t_0,t_{i}))
% CONVENTIONS:
%   - delta(.,.) : 30/360
%   - TO BE COMPLETED...

%% BOOTSTRAP PROCEDURE ----------------------------------------------------
% ----------------------- First iteration (i = 1): ------------------------
% Compute the initial survival probability.
% The formula neglects accrual term and adjusts only for recovery & spread.
survProbs(1) = (1 - recovery) /...
               (1 - recovery + spreadsCDS(1) * delta_t(1));

% Calculate the first hazard rate from the initial survival probability.
intensities(1) = -log(survProbs(1)) / delta_t(1);

% ----------------------- Next iterations (i > 1): ------------------------
% ---------------------- For subsequent maturities ------------------------

% Initialize the cumulative variables:
% RPV01 accumulates the discounted survival probabilities.
RPV01 = delta_t(1) * survProbs(1) * discounts_CDS(1);
% e accumulates the discounted expected loss (1 - survival probability).
e = discounts_CDS(1) * (1 - survProbs(1));

% Loop over each i >= 2 CDS maturity to bootstrap survival probabilities.
for i = 2:N
    % Compute the NUMERATOR of the survival probability formula:
    % It combines the discounted survival from the previous period with 
    % adjustments for EL and RPV01 weighted by the CDS spread.
    NUM = (1 - recovery) * discounts_CDS(i) * survProbs(i-1) +...
          (1 - recovery) * e - spreadsCDS(i) * RPV01;
    % e = SUM_{j=1}^{i-1} B(t_0,t_j) [P(t_0,t_{j-1}) - P(t_0,t_j)]
    % (cf notations above in the theory explanation)
    
    % Compute the DENOMINATOR of the formula:
    % This adjusts the discounted survival probability for the spread over
    % the time increment.
    DEN = ((1 - recovery) + spreadsCDS(i) * delta_t(i)) * discounts_CDS(i);
    
    % Calculate the survival probability for the current maturity.
    survProbs(i) = NUM / DEN;
    
    % Update RPV01 by adding the current period's contribution, product of:
    % time increment, new survival probability, & corresponding DF.
    RPV01 = RPV01 + delta_t(i) * survProbs(i) * discounts_CDS(i);
    
    % Update the expected loss accumulator 'e' by adding the discounted 
    % change in survival probability.
    % Remember that:
    % e = SUM_{j=1}^{i-1} B(t_0,t_j) [P(t_0,t_{j-1}) - P(t_0,t_j)]
    % cf notations above in the theory explanation.
    e = e + discounts_CDS(i) * (survProbs(i-1) - survProbs(i));
    
    % Compute the hazard intensity for the interval [i - 1, i]:
    % cf formula in the theory explanation above.
    intensities(i) = log(survProbs(i-1) / survProbs(i)) / delta_t(i);
end

end