function varargout = parcorr(Series , nLags , P , nSTDs)
%PARCORR Compute or plot sample partial auto-correlation function.
%   Compute or plot the sample partial auto-correlation function (partial ACF) 
%   of a univariate, stochastic time series. The partial ACF is computed by 
%   fitting successive autoregressive models of orders 1,2, ... by ordinary 
%   least squares, retaining the last coefficient of each regression. When 
%   called with no output arguments, PARCORR displays the sequence with 
%   confidence bounds.
%
%   [PartialACF, Lags, Bounds] = parcorr(Series)
%   [PartialACF, Lags, Bounds] = parcorr(Series , nLags , R , nSTDs)
% 
%   Optional Inputs: nLags , R , nSTDs
%
% Inputs:
%  Series - Vector of observations of a univariate time series for which the
%     sample partial ACF is returned or plotted. The last row of Series 
%     contains the most recent observation of the stochastic sequence.
%
% Optional Inputs:
%   nLags - Positive, scalar integer indicating the number of lags of the 
%     partial ACF to compute. If empty or missing, the default is to compute 
%     the partial ACF sequence at lags 0,1,2, ... T, where T is equal to the 
%     minimum[20 , length(Series)-1].
%
%   R - Non-negative integer scalar indicating the number of lags beyond which 
%     the theoretical partial ACF is assumed zero. Under the hypothesis that 
%     the underlying Series is really an AR(R) process, the estimated partial 
%     ACF coefficients at lags > R are approximately zero-mean, independently
%     distributed Gaussian variates. In this case, the standard error of the 
%     estimated partial ACF coefficients of a fitted Series with N observations
%     is approximately 1/sqrt(N) for lags > R. If R is empty or missing, the 
%     default is R = 0. R must be less than nLags.
%
%   nSTDs - Positive scalar indicating the number of standard deviations of the 
%     sample partial ACF estimation error to display assuming that Series is
%     an AR(R) process. If the Rth regression coefficient (i.e., the last OLS 
%     regression coefficient of Series regressed on a constant and R of its 
%     lags) is fitted with N observations, specifying nSTDs will result in 
%     confidence bounds at +/-(nSTDs/sqrt(N)). If empty or missing, default is
%     nSTDs = 2 (i.e., approximate 95% confidence interval).
%
% Outputs:
%   PartialACF - Sample partial ACF of Series. PartialACF is a vector of length
%     nLags + 1 corresponding to lags 0,1,2,...,nLags. The first element of 
%     PartialACF is defined to be unity (i.e., PartialACF(1) = 1 = OLS 
%     regression coefficient of Series regressed upon itself), and is included 
%     as a reference.
%
%   Lags - Vector of lags corresponding to PartialACF (0,1,2,...,nLags).
%
%   Bounds - Two element vector indicating the approximate upper and lower
%     confidence bounds assuming that Series is an AR(R) process. Note that 
%     Bounds is approximate for lags > R only.
%
% Example:
%   Create a stationary AR(2) process from a sequence of 1000 Gaussian deviates,
%   then visually assess whether the partial ACF is zero for lags > 2:
%
%     randn('state',0)               % Start from a known state.
%     x = randn(1000,1);             % 1000 Gaussian deviates ~ N(0,1).
%     y = filter(1,[1 -0.6 0.08],x); % Create a stationary AR(2) process.
%     parcorr(y , [] , 2)            % Inspect the P-ACF with 95% confidence.
%
% See also CROSSCORR, AUTOCORR, FILTER.

%   Copyright 1999-2003 The MathWorks, Inc.   
%   $Revision: 1.6.2.2 $  $Date: 2007/09/11 11:46:47 $

%
% References:
%   Box, G.E.P., Jenkins, G.M., Reinsel, G.C., "Time Series Analysis: 
%     Forecasting and Control", 3rd edition, Prentice Hall, 1994.
%   Hamilton, J.D., "Time Series Analysis", Princeton University Press, 1994.
%

%
% Ensure the sample data is a VECTOR.
%

[rows , columns]  =  size(Series);

if (rows ~= 1) && (columns ~= 1) 
    error('econ:parcorr:NonVectorInput' , ' Input ''Series'' must be a vector.');
end

rowSeries   =  size(Series,1) == 1;

Series      =  Series(:);       % Ensure a column vector
n           =  length(Series);  % Raw sample size.
defaultLags =  20;              % BJR recommend about 20 lags for partial ACFs.

%
% Ensure the number of lags, nLags, is a positive 
% integer scalar and set default if necessary.
%

if (nargin >= 2) && ~isempty(nLags)
   if numel(nLags) > 1
      error('econ:parcorr:NonScalarLags' , ' Number of lags ''nLags'' must be a scalar.');
   end
   if (round(nLags) ~= nLags) || (nLags <= 0)
      error('econ:parcorr:NonPositiveIntegerLags' , ' Number of lags ''nLags'' must be a positive integer.');
   end
   if nLags > (n - 1)
      error('econ:parcorr:LagsTooLarge' , ' Number of lags ''nLags'' must not exceed ''Series'' length - 1.');
   end
else
   nLags  =  min(defaultLags , n - 1);
end

%
%  Ensure the hypothesized number of lags, P, is a non-negative integer
%  scalar, and set default if necessary.
%

if (nargin >= 3) && ~isempty(P)
   if numel(P) > 1
      error('econ:parcorr:NonScalarP' , ' Number of lags ''P'' must be a scalar.');
   end
   if (round(P) ~= P) || (P < 0)
      error('econ:parcorr:NegativeIntegerP' , ' Number of lags ''P'' must be a non-negative integer.');
   end
   if P >= nLags
      error('econ:parcorr:PTooLarge' , ' ''P'' must be less than ''nLags''.');
   end
else
   P  =  0;       % Set default.
end

%
%  Ensure the number of standard deviations, nSTDs, is a positive 
%  scalar and set default if necessary.
%

if (nargin >= 4) && ~isempty(nSTDs)
   if numel(nSTDs) > 1
      error('econ:parcorr:NonScalarSTDs' , ' Number of standard deviations ''nSTDs'' must be a scalar.');
   end
   if nSTDs < 0
      error('econ:parcorr:NegativeSTDs' , ' Number of standard deviations ''nSTDs'' must be non-negative.');
   end
else
   nSTDs =  2;     % Default is 2 standard errors (~95% condfidence interval).
end

%
% Create a lagged regression matrix & allocate storage for the partial ACF.
%

X          =  lagmatrix(Series , [1:nLags]);
partialACF =  [1 ; zeros(nLags , 1)];

%
% Compute partial ACF by fitting successive order AR models 
% by OLS, retaining the last coefficient of each regression.
%

for order = 1:nLags
   [Q , R]             =  qr([ones((length(Series)-order),1)  X(order+1:end,1:order)] , 0);
   b                   =  R\(Q'*Series(order+1:end));
   partialACF(order+1) =  b(end);
end

%
% Compute approximate confidence bounds using the Box-Jenkins-Reinsel 
% approach, equations 3.2.36 and 6.2.3, on pages 68 and 188, respectively. 
%
% Note a subtle point here: The Pth autoregressive model 'fit' via OLS 
% makes use of only the most recent (n - P) observations. Since the 
% approximate confidence bounds for the hypothesized P is of interest 
% only for lags > P, and the (P+1)th AR model uses (n - (P + 1) = n - p - 1
% observations, the 'n' in BJR equation 3.2.36 (i.e., the number of 
% observations used in 'fitting') is taken to be (n - P - 1) rather than
% the original length of Series. Moreover, the effective number of 
% observations used in 'fitting' each successive AR model will decrease 
% by one observation for each lag. For even moderate sample sizes, this
% approximation should make little difference.
%

bounds  =  [nSTDs ; -nSTDs] ./ sqrt(n - P - 1);
Lags    =  [0:nLags]';

if nargout == 0

%
%  Plot the sample partial ACF. Note the partial ACF at lag 0 is defined to be 1.
%
   lineHandles  =  stem(Lags , partialACF , 'filled' , 'r-o');
   set   (lineHandles(1) , 'MarkerSize' , 4)
   grid  ('on')
   xlabel('Lag')
   ylabel('Sample Partial Autocorrelations')
   title ('Sample Partial Autocorrelation Function')
   hold  ('on')
%
%  Plot the confidence bounds under the hypothesis that the underlying 
%  Series is really an AR(P) process. The following approximation gives
%  an indication of whether the partial ACF is effectively zero beyond 
%  lag P. For this reason, the confidence bounds (horizontal lines) appear 
%  over the partial ACF ONLY for lags > P (i.e., P+1, P+2, ... nLags).
%  In other words, the confidence bounds enclose ONLY those lags for 
%  which the null hypothesis is assumed to hold.
%

   plot([P+0.5 P+0.5 ; nLags nLags] , [bounds([1 1]) bounds([2 2])] , '-b');
   plot([0 nLags] , [0 0] , '-k');
   hold('off')

   if max(partialACF) <= 1
      a  =  axis;
      axis([a(1:3) 1]);
   end

else

%
%  Re-format outputs for compatibility with the SERIES input. When SERIES is
%  input as a row vector, then pass the outputs as a row vectors; when SERIES
%  is a column vector, then pass the outputs as a column vectors.
%
   if rowSeries
      partialACF =  partialACF.';
      Lags       =  Lags.';
      bounds     =  bounds.';
   end

   varargout  =  {partialACF , Lags , bounds};

end
