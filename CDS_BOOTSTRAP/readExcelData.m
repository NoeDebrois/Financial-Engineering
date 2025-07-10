function [dates, rates] = readExcelData(filename, formatData)
    % Reads data from Excel
    %  It reads bid/ask prices and relevant dates
    %  All input rates are in % units
    %
    % INPUTS:
    %  filename: Excel file name where data are stored
    %  formatData: data format in Excel
    % 
    % OUTPUTS:
    %  dates: struct with settlementDate, deposDates, futuresDates, swapDates
    %  rates: struct with deposRates, futuresRates, swapRates
    
    %% Dates from Excel
    
    % Settlement date
    settlement = readtable(filename, 'Range', 'E8:E8', 'ReadVariableNames', false);
    % Date conversion
    dates.settlement = datenum(settlement{1,1});
    
    % Dates relative to deposits
    date_depositi = readtable(filename, 'Range', 'D11:D16', 'ReadVariableNames', false);
    dates.depos = datenum(date_depositi.Var1);
    
    % Dates relative to futures: calc start & end
    date_futures_read = readtable(filename, 'Range', 'Q12:R20', 'ReadVariableNames', false);
    numberFutures = size(date_futures_read, 1);
    
    dates.futures = ones(numberFutures, 2);
    dates.futures(:, 1) = datenum(date_futures_read.Var1);
    dates.futures(:, 2) = datenum(date_futures_read.Var2);
    
    % Dates relative to swaps: expiry dates
    date_swaps = readtable(filename, 'Range', 'D39:D55', 'ReadVariableNames', false);
    dates.swaps = datenum(date_swaps.Var1);
    
    %% Rates from Excel (Bids & Asks)
    
    % Depos
    tassi_depositi = readtable(filename, 'Range', 'E11:F16', 'ReadVariableNames', false);
    rates.depos = tassi_depositi{:,:} / 100;
    
    % Futures
    tassi_futures = readtable(filename, 'Range', 'E28:F36', 'ReadVariableNames', false);
    % Rates from futures
    tassi_futures = 100 - tassi_futures{:,:};
    rates.futures = tassi_futures / 100;
    
    % Swaps
    tassi_swaps = readtable(filename, 'Range', 'E39:F55', 'ReadVariableNames', false);
    rates.swaps = tassi_swaps{:,:} / 100;
    
    end % readExcelData
    