% This is a function that takes two corresponding time series and perofmrs
% an intense validation following the proceedure described by Chris A.
% Gueymard in his 2014 paper:
%
% Christian A. Gueymard. 2014. A review of validation methodologies and 
% statistical performance indicators for modeled solar radiation data: 
% Towards a better bankability of solar projects. Journal of Renewable and 
% Sustainable Energy Reviews. Volume 39. Pages 1024-1034.
%
% -----------------------------------------------------------------------
%                               CONTENTS
% -----------------------------------------------------------------------
% The method is categorised into four classes labeled A, B C and D.
% 
% A. Class A—indicators of dispersion
%  A.1. Mean bias error (MBE)
%  A.2. Root mean square error (RMSE). 
%  A.3. Mean absolute error (MAE)
%  A.4. Standard deviation of the residual (SD)
%  A.A. Coefficient of determination (R2)
%  A.6. Slope of best-fit line (SBF)
%  A.7. Uncertainty at 95% (U95)
%  A.8. t-statistic (TS)
%
% B Class B—indicators of overall performance.
%  B.1. Nash–Sutcliffe's efficiency (NSE)
%  B.2. Willmotts's index of agreement (WIA)
%  B.3. Legates's coefficient of efficiency (LCE)
%
% C. Class C—indicators of distribution similitude
%  C.1. Kolmogorov-Smirnov test integral 
%  C.2. OVER statistic of relative frequency of exeedence situations (OVER)
%
% D. Class D—visual indicators
%  D.1. Taylor diagram
%  D.2. Mutual information diagram
%  D.3. Boxplot Diagram
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% Observations and Predictions must be numerical and of the same size.
% Observations is treated as the correct value. They should be column
% vectors.
%
% class_selection must be a cell indicating which Class tests are to be
% performed. The default tests are all A, B, C and D. To indicate an
% individual test, set class_selection to {'B'} for example.
%
% -----------------------------------------------------------------------
%                               OUTPUTS
%  -----------------------------------------------------------------------
% The output is a struct containing the results, depending on the selected 
% class_selection. Within each struct contains the respective results. 
% For example, validation_struct.MBE contains the mean bias error of the
% validation.



function validation_struct = FourClassValidation(Observations,Predictions,class_selection)
 %% preliminary checks and preprocessing
% extract to simple variable names
O = Observations;
P = Predictions;

% default class selection to all
if ~exist('class_selection','var')
    class_selection = {'A','B','C','D'};
end


% checks
if (~isnumeric(O) || ~isnumeric(P))
    error('Observations and Predictions must be of class numeric')
end
if size(O)~=size(P)
    error('Observations and Predictions must be of equal size')
end
if (min(size(O))~=1 || min(size(P))~=1)
    error('Observations and Predictions must be a single corresponding vector')
end

% insist on column vectors
O = reshape(O,[numel(O),1]);
P = reshape(P,[numel(P),1]);

% removal of NaN values
not_nan_inds = (~isnan(O) & ~isnan(P));
O = O(not_nan_inds);
P = P(not_nan_inds);

% define common usages
 Om = mean(O);
 Pm = mean(P);
 N = length(O);
 
 validation_struct.Om = Om;
 validation_struct.Pm = Pm;
 validation_struct.N = N;

%% Class A - indicators of dispersion
% These are the indicators that the majority of readers should be most 
% familiar with. They are all expressed here in percent (of Om) rather than
% in absolute units (W/m2 for irradiances, or MJ/m2 or kWh/ m2 for 
% irradiations) because non-expert stakeholders can much more easily 
% understand percent results. In any case, stating the value of Om in all 
% validation results allows the experts to convert back the percent figures 
% into absolute units if they so desire. Formulas in this section are well 
% established and do not need further references.

if max(strcmpi(class_selection,'A')==1)
%-------------------------------------------------------------------------    
    % A.1  Mean bias error (MBE)
    validation_struct.MBE = 100./Om .* sum(P-O);
%-------------------------------------------------------------------------    
    % A.2 Root mean square error (RMSE)
    validation_struct.RMSE = 100./Om .* ( sum(P-O).^2 ./ N ).^0.5;
%-------------------------------------------------------------------------    
    % A.3 Mean absolute error (MAE)
    validation_struct.MAE = 100./Om .* sum(abs(P-O));
%-------------------------------------------------------------------------    
    % A.4 Standard deviation of the residual (SD)
    validation_struct.SD = 100./Om .* ( sum(N.*(P-O).^2)  - sum(P-O).^2 ).^0.5 ./ N;
%-------------------------------------------------------------------------    
    % A.5 Coefficient of determination (R2)
    validation_struct.R2 = (  sum((P-Pm) .* (O-Om)) ./ sum((P-Pm).^2 .* (O-Om).^2)  ).^2;
%-------------------------------------------------------------------------    
    % A.6 Slope of best-fit line (SBF)
    validation_struct.SBF = sum((P-Pm).*(O-Om)) ./ sum(O-Om).^2;
%-------------------------------------------------------------------------    
    % A.7 Uncertainty at 95%
    validation_struct.U95 = 1.96 .* (validation_struct.SD.^2 + validation_struct.RMSE.^2).^0.5;
%-------------------------------------------------------------------------    
    % A.8 t-statistic (TS)
    validation_struct.TS = ( (N-1) .* validation_struct.MBE.^2 ./ (validation_struct.RMSE.^2-validation_struct.MBE.^2) ).^0.5;    
%-------------------------------------------------------------------------
end


%% Class B - Indicators of overall performance
% These are indicators that are less common in the solar field than those
% of Class A. They convey relatively similar information as those of Class 
% A, with the cosmetic advantage that a higher value indicates a better
% model.

if max(strcmpi(class_selection,'B')==1)
%-------------------------------------------------------------------------
    % B.1 Nash-Sutcliffe's efficiency (NSE)
    validation_struct.NSE = 1 - sum(P-O).^2 ./ sum(O-Om).^2;
%-------------------------------------------------------------------------
    % B.2 Willmotts's index of agreement (WIA)
    validation_struct.WIA = 1 - sum(P-O).^2 ./ sum(abs(P-Om) + abs(O-Om)).^2;
%-------------------------------------------------------------------------    
    % B.3 Lagates's coefficient of efficiency (LCE)
    validation_struct.LCE = 1 - sum(abs(P-O)) ./ sum(abs(O-Om));
%-------------------------------------------------------------------------    
% LSE and NSE vary between 1 for perfect agreement and -inf for complete
% disagreement, whereas WIE varies only between 1 and 0.        

end

%% Class C
% The goal is to compare one or more cumulative frequency distribution
% of modeled data to that of a reference dataset. Can one or more single 
% number provide a measure of the similitude between two or more 
% distributions? Substantial progress in that direction resulted from an
% initial study by Polo et al. [1],who proposed to use the Kolmogorov–
% Smirnov test when comparing different cumulative distribution functions 
% (CDFs), because of its advantage of being non- parametric and valid for 
% any kind of CDF. Espinar et al. [2] developed the method further, now
% referring to it as the Kolmogorov–Smirnov test Integral (KSI)

if max(strcmpi(class_selection,'C')==1)
%-------------------------------------------------------------------------   
    % C.1 Kolmogorov-Smirnov test Integral (KSI)
    % irradiance must be binned into x by intervals of n
    n = 100;
    xbins = (0:n:1500);
    xmin = min(xbins);
    xmax = max(xbins);
    Od = histc(O,xbins);
    Pd = histc(P,xbins);            
    % absolute difference between the two normalised distributions
    Dn = abs(Od - Pd); 
    % pure function of N obtained from [3], though simplified to constant \approx 1.63.
    Dc = 1.63; 
    Ac = Dc .* (xmax - xmin);
    fun = @(x) Dn(x);
    validation_struct.KSI = 100/Ac .* integral(fun,xmin,xmax);
    % KSI is 0 if the two distributions being compared can be considered
    % identical in a statistical sense.
%-------------------------------------------------------------------------    
    % C.2 Relative frequency of exeedence situations (OVER)
    % when the normalised distribution of modelled data points in specific
    % bins exceeds the critical limit that would make it statistically
    % undistinguisable from the reference distribution.
    fun = @(x) max([Dn(x)-Dc,0]);
    validation_struct.OVER = 100/Ac .* integral(fun,xmin,xmax);
    % OVER is 0 if the normalised distribution always remains below Dc.
    % OVER can be null indicating that the distribution of the  predictions
    % generally respect those of the predictions.
    
%-------------------------------------------------------------------------    
    % C.3 Combined Performance Index (CPI)
    % The interest of CPI is that it combines conventional information
    % about dispersion and biase (through RMSE) with information about
    % distribution likenesses (through KSI and OVER), whilst maintaining a
    % high degree of discrimination between the different models. This
    % feature is of course highly desireable when comparing differnet
    % models of similar performance. This is arguably the most significant
    % statistic to compare different model performance. 
    
    % The CPI requires the RMSE of the two indices, should class_section A
    % not have been performed, the RMSE must also be calculated
    if max(strcmpi(class_selection,'A'))==0
          % A.2 Root mean square error (RMSE)
          validation_struct.RMSE = 100./Om .* ( sum(P-O).^2 ./ N ).^0.5;
    end
    % all values must be in percentages.
    validation_struct.CPI = (validation_struct.KSI + OVER + 2.*validation_struct.RMSE)./4;
%-------------------------------------------------------------------------
end

%% Class D
% This category is completely different from the three previous ones 
% because the goal here is to obtain a visualization rather than summary 
% statistics in the form of a few numbers
if max(strcmpi(class_selection,'D')==1)
     % The first recommended plots for class D is a Taylor diagram detailed 
     % by KE Taylor [4] that combines RMSE, SD and R2 into a single polar
     % diagram. It is ideal for comparing the performance of many different
     % models. 
     
     % The second suggestion is a Mutual Information Diagram, which is a
     % revision of the Taylor diagram proposed by [5].
     
     % The box plot is another decent variation for demonstrating
     % performance at different sites.
    
    
end

end

%% References
%
% [1] Polo J, Zarzalejo LF, Ramirez L, Espinar B. Iterative filtering of 
% ground data for qualifying statistical models for solar irradiance 
% estimation from satellite data. Sol. Energy 2006; 80:240–7
%
% [2] Espinar B, Ramirez L, Drews A, Beyer HG, Zarzalejo LF, Polo J, et 
% al. Analysis of different comparison parameters applied to solar 
% radiation data from satellite and German radiometric stations. Sol Energy
% 2009;83:118–25.
%
% [3] Marsaglia G, Tsang WW, Wang J. Evaluating Kolmogorov's Distribution.
% J Stat Softw 2003;8:1–4
%
% [4] Taylor KE. Summarizing multiple aspects of model performance in a 
% single diagram. J Geophys Res 2001;106D:7183–92.
%
% [5] Correa CD, Lindstrom P. The mutual information diagram for 
% uncertainty visualization. Int J Uncertain Quantif 2013;3:187–201.
