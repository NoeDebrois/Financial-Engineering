function [datesCDS, survProbs, intensities] = bootstrapCDS_accrual(datesDF, discounts, datesCDS, spreadsCDS, recovery)
% bootstrapCDS_accrual computes the survival probabilities and hazard rates
% (intensities) for a set of CDS instruments, including an accrual term.
%
% INPUT:
%   datesDF   : Vector of dates corresponding to discount factors
%   discounts : Vector of discount factors
%   datesCDS  : Vector of CDS maturity/expiry dates
%   spreadsCDS: Vector of CDS spreads
%   recovery  : Recovery rate for the reference entity
%
% OUTPUT:
%   datesCDS   : Same as input, dates of expiry of the CDS instruments
%   survProbs  : Survival probabilities at each CDS date
%   intensities: Hazard rates (intensities) at each CDS date
%
% DESCRIPTION:
%   This function bootstraps the survival probabilities and hazard rates
%   using a piecewise-constant hazard rate model, explicitly including the
%   accrual on the premium leg of the CDS. The accrual term accounts for 
%   the fraction of the premium that would be paid if default occurs 
%   between payment dates.

%% PRELIMINARY COMPUTATIONS -----------------------------------------------
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
                    /100 .* dates_CDS_zrates );

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
% [E] is the same as explained in bootstrapCDS_NOaccrual.m except for an
% added term on the premium leg:
% [E]
% S_{iY} * SUM_{j=1}^{i} delta(t_{j-1},t_j) B(t_0,t_j) P(t_0,t_j)
% +
% S_{iY} * SUM_{j=1}^{i} delta(t_{j-1},t_j)/2 B(t_0,t_j) [P(t_0,t_{j-1}) -
%                                                         P(t_0,t_j)]
% = 
% (1 - RR) * SUM_{j=1}^{i} B(t_0,t_j) [P(t_0,t_{j-1}) - P(t_0,t_j)]
% i.e., PREMIUM LEG for maturity i Y = PROTECTION LEG for maturity i Y.
% NB: take a look at explanation in bootstrapCDS_NOaccrual.m if needed.
% NB2: for i>1, rewrite [E] the following way:
%      P(t_0,t_i) = NUM / DEN where:
%      NUM = [1 - RR - delta(t_{i-1},t_i)/2 * S_{iY}] * B(t_0,t_i) *
%            P(t_0,t_{i-1}) + [1 - RR] * e(t_{i-1})
%            - S_{iY} * [RPV01(t_{i_1}) + e_cum_deltas(t_{i-1})]
%      DEN = [1 - recovery + delta(t_{i-1},t_i)/2 * S_{iY}] * B(t_0,t_i)
%      where:
%      e(t_{i-1} = SUM_{j=1}^{i-1} B(t_0,t_j) [P(t_0,t_{j-1}) - P(t_0,t_j)]
%      e_cum_deltas(t_{i-1}) = SUM_{j=1}^{i-1} delta(t_{j-1},t_j)/2 * 
%                              B(t_0,t_j) * [P(t_0,t_{j-1}) - P(t_0,t_j)]
%      RPV01(t_{i_1}) = SUM_{j=1}^{i-1} delta(t_{j-1},t_j) * B(t_0,t_j) * 
%                       P(t_0,t_j)
% CONVENTIONS:
%   - delta(.,.) : 30/360
%   - TO BE COMPLETED...

%% BOOTSTRAP PROCEDURE ----------------------------------------------------
% ----------------------- First iteration (i = 1): ------------------------
% The formula for the first survival probability includes an accrual term 
% proportional to (delta_t(1)/2). The approach is an approximation where 
% the default is assumed to happen in the middle of the period on average.

survProbs(1) = (1 - recovery - spreadsCDS(1) * delta_t(1)/2) / ...
               (1 - recovery + spreadsCDS(1) * delta_t(1)/2);

% Calculate the first hazard rate from the initial survival probability.
intensities(1) = -log(survProbs(1)) / delta_t(1);

% ----------------------- Next iterations (i > 1): ------------------------
% ---------------------- For subsequent maturities ------------------------

% Initialize the cumulative variables:
% RPV01 accumulates the discounted survival probabilities.
RPV01 = delta_t(1) * survProbs(1) * discounts_CDS(1);
% e accumulates the discounted expected loss (1 - survival probability).
e = discounts_CDS(1) * (1 - survProbs(1));
% e_cum_deltas tracks the discounted accrual portion.
e_cum_deltas = delta_t(1)/2 * discounts_CDS(1) * (1 - survProbs(1));

% Loop over each i >= 2 CDS maturity to bootstrap survival probabilities.
for i = 2:N
    % Compute the NUMERATOR of the survival probability formula:
    % NUM = [1 - RR - delta(t_{i-1},t_i)/2 * S_{iY}] * B(t_0,t_i) *
    %       P(t_0,t_{i-1}) + [1 - RR] * e(t_{i-1})
    %       - S_{iY} * [RPV01(t_{i_1}) + e_cum_deltas(t_{i-1})]
    % (cf notations above in the theory explanation)
    NUM = (1 - recovery) * e - spreadsCDS(i) * RPV01... 
          - spreadsCDS(i) * e_cum_deltas +...
          (1 - recovery - spreadsCDS(i) * delta_t(i)/2)...
          * discounts_CDS(i) * survProbs(i-1);

    % Compute the DENOMINATOR of the formula:
    % DEN = [1 - recovery + delta(t_{i-1},t_i)/2 * S_{iY}] * B(t_0,t_i)
    % (cf notations above in the theory explanation)
    DEN = (1 - recovery + spreadsCDS(i) * delta_t(i)/2) * discounts_CDS(i);
    
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
    
    % Update e_cum_deltas: accrual part
    % Remember that:
    % e_cum_deltas = SUM_{j=1}^{i-1} delta(t_{j-1},t_j)/2 B(t_0,t_j) 
    %                [P(t_0,t_{j-1}) - P(t_0,t_j)]
    e_cum_deltas = e_cum_deltas + (delta_t(i)/2) * discounts_CDS(i) * ...
                   (survProbs(i-1) - survProbs(i));
    
    % Compute the hazard intensity for the interval [i - 1, i]:
    % cf formula in the theory explanation above.
    intensities(i) = log(survProbs(i-1) / survProbs(i)) / delta_t(i);
end

end