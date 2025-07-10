%% BOOTSTRAP OF THE DISCOUNT CURVE ----------------------------------------
formatData = 'dd/mm/yyyy'; % Date format string.

% Reads market data from an Excel file:
% Note: this function has been modified to work in MacOS.
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);

% The bootstrap fct returns a set of dates (including settlement) & the 
% corresponding DFs computed from mkt data:
[dates, discounts] = bootstrap(datesSet, ratesSet);
z_Rates = zeroRates(dates, discounts); % DF -> Zero Rates.
% END OF BOOTSTRAP OF THE DISCOUNT CURVE ----------------------------------





%% ASSIGNEMENT 3 - EXERCISE 2 - CDS Bootstrap -----------------------------

% Input data for the CDS: ----------------- [A MODIFIER EN FCT DE L'ENONCE]
recovery = 40/100; % Recovery rate.

spreadsCDS = [40 44 47 49 51 52]' * 1e-4; % CDS spreads in decimal for 
% maturities in datesSet.swaps.

datesCDS = datesSet.swaps(1:6); % CDS maturities from mkt data structure.
% -------------------------------------------------------------------------

% [1] CDS BOOTSTRAPPING NEGLECTING ACCRUAL TERM ---------------------------
% [APPROXIMATION METHOD -> flag = 1]
% Lambda is assumed piecewise constant.
flag = 1;
[datesCDS, survProbs_Noaccural, intensities_Noaccrual] =... 
    bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, flag, recovery);
% bootstrapCDS.m -> bootstrapCDS_NOaccrual.m
% bootstrapCDS_NOaccrual.m is fully commented.

fprintf('Survival Probabilities, no accrual:\n');
fprintf('%f\n', survProbs_Noaccural);
fprintf('Intensities, no accrual:\n');
fprintf('%f\n', intensities_Noaccrual);
% [1] END OF CDS BOOTSTRAPPING NEGLECTING ACCRUAL TERM --------------------


% [2] CDS BOOTSTRAPPING INCLUDING ACCRUAL TERM ----------------------------
% [EXACT METHOD -> flag = 2]
% Lambda is assumed piecewise constant.
flag = 2;
[datesCDS, survProbs_accrual, intensities_accrual] =...
    bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, flag, recovery);
% bootstrapCDS.m -> bootstrapCDS_accrual.m

fprintf('Survival Probabilities, with accrual:\n');
fprintf('%f\n', survProbs_accrual);
fprintf('Intensities, with accrual:\n');
fprintf('%f\n', intensities_accrual);
% [2] END OF CDS BOOTSTRAPPING INCLUDING ACCRUAL TERM ---------------------


% COMPARISON BETWEEN METHOD 1 & METHOD 2:
fprintf("Max surv probability error neglecting accrual term: %.2d \n",... 
        max(abs(survProbs_accrual - survProbs_Noaccural)));
fprintf("Max intensity error neglecting accrual term: %.2d \n",...
        max(abs(intensities_accrual - intensities_Noaccrual)));
% The accrual term can be in general neglected, as it is 2 orders of 
% magnitude smaller than the other terms.


% [3] CDS BOOTSTRAPPING USING JT APPROACH ---------------------------------
% [JT-APPROXIMATION METHOD -> flag = 3]
% Lambda is assumed constant.
flag = 3;
[datesCDS, survProbs_JT, intensities_JT] =...
    bootstrapCDS(dates, discounts, datesCDS, spreadsCDS, flag, recovery);
% bootstrapCDS.m -> bootstrapCDS_JT.m

fprintf('Survival Probabilities, JT approximation:\n');
fprintf('%f\n', survProbs_JT);
fprintf('Intensities, JT approximation:\n');
fprintf('%f\n', intensities_JT);
% [3] END OF CDS BOOTSTRAPPING USING JT APPROACH --------------------------

% COMPARISON BETWEEN ...
fprintf("\nMax prob error with JT approximation: %.2d \n",...
        max(abs(survProbs_accrual - survProbs_JT)));
fprintf("Max intensities error with JT approximation: %.2d \n\n",...
        max(abs(intensities_accrual - intensities_JT)));
% [3] END OF CDS BOOTSTRAPPING USING JT APPROACH --------------------------


% PLOTS -------------------------------------------------------------------
% Plot survival probabilities for the three methods
figure;
x = (0:6);  % Time index: 0 represents time 0 (settlement), then each year
plot(x, [1; survProbs_Noaccural], '-s', 'LineWidth', 1.5, 'MarkerSize',10);
hold on;
plot(x, [1; survProbs_accrual], '-o', 'LineWidth', 1.1, 'MarkerSize', 5);
plot(x, [1; survProbs_JT], '-^', 'LineWidth', 1.1, 'MarkerSize', 5);
legend('Approximate', 'Exact', 'JT', 'Location', 'northeast');
set(gca, 'YGrid', 'on');
title('Survival Probabilities', 'FontSize', 14);
xlabel('Time', 'FontSize', 12);
ylabel('Probability', 'FontSize', 12);
xticklabels({'02/02/23','1y','2y','3y','4y','5y','6y'});
xtickangle(45);
hold off;

% Plot intensities for the approximate and exact methods
figure;
x_i = linspace(0, 6, 500);
% Interpolate intensities for a smoother curve; using "next" interpolation 
% method:
plot(x_i,...
     interp1(x, [intensities_Noaccrual(1); intensities_Noaccrual], x_i,...
     "next", "extrap"),...
     '-s', 'LineWidth', 1.5, 'MarkerSize', 0.75);
hold on;
plot(x_i,...
     interp1(x, [intensities_accrual(1); intensities_accrual], x_i,...
     "next", "extrap"),...
     '-o', 'LineWidth', 1.1, 'MarkerSize', 0.55);
plot(x_i,... 
    interp1(x, [intensities_JT(1); intensities_JT], x_i,...
    "next","extrap"),...
    '-^', 'LineWidth', 1.1, 'MarkerSize', 0.55);
legend('Approximate', 'Exact','JT ','Location', 'southeast');
set(gca, 'YGrid', 'on');
title('Intensities', 'FontSize', 14);
xlabel('Time', 'FontSize', 12);
ylabel('Intensity', 'FontSize', 12);
xticklabels({'02/02/23','1y','2y','3y','4y','5y','6y'});
xtickangle(45);
hold off;
% -------------------------------------------------------------------------