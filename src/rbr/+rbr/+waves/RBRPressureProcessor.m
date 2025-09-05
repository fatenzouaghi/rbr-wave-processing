@Date    : 2025-09-02
@Author  : Faten Zouaghi (faten_zouaghi@uqar.ca)

classdef RBRPressureProcessor
% RBRPressureProcessor
% End-to-end pipeline: RBR pressure (.rsk) + meteo CSV -> wave time series.

properties (Constant)
    rho   = 1023;    % seawater density [kg/m^3]
    g     = 9.81;    % gravity [m/s^2]
    R_gas = 8.314;   % molar gas constant [J/(mol*K)]
    M_air = 0.02896; % dry-air molar mass [kg/mol]
end

% --- No site-specific defaults: Specify values per deployment/data.
properties
    fs        (1,1) double = NaN;   % sampling frequency [Hz]
    nfft      (1,1) double = NaN;   % FFT length
    minFreq   (1,1) double = NaN;   % Hz
    maxFreq   (1,1) double = NaN;   % Hz
    igCutoff  (1,1) double = NaN;   % Hz (IG vs sea–swell split)
    windowSec (1,1) double = NaN;   % seconds (block length)
    critWL    (1,1) double = NaN;   % m (min WL to accept a spectrum)
end

methods
    function out = run(obj, varargin)
        ip = inputParser;
        addParameter(ip,'RSK','',@(s)ischar(s)||isstring(s));
        addParameter(ip,'MeteoCSV','',@(s)ischar(s)||isstring(s));
        addParameter(ip,'alti',[],@isnumeric);
        addParameter(ip,'zmembrane',[],@isnumeric);
        addParameter(ip,'zbottom',[],@isnumeric);
        addParameter(ip,'fs',[],@isnumeric);
        addParameter(ip,'nfft',[],@isnumeric);
        addParameter(ip,'minFreq',[],@isnumeric);
        addParameter(ip,'maxFreq',[],@isnumeric);
        addParameter(ip,'igCutoff',[],@isnumeric);
        addParameter(ip,'windowSec',[],@isnumeric);
        addParameter(ip,'critWL',[],@isnumeric);
        parse(ip,varargin{:});
        A = ip.Results;

  
        if isempty(A.RSK),      A.RSK      = getenv('RBR_RSK');    end
        if isempty(A.MeteoCSV), A.MeteoCSV = getenv('RBR_METEO');  end
        assert(~isempty(A.RSK)      && exist(A.RSK,'file')==2, 'RSK file not found.');
        assert(~isempty(A.MeteoCSV) && exist(A.MeteoCSV,'file')==2, 'Meteo CSV not found.');
        assert(~isempty(A.alti) && ~isempty(A.zmembrane) && ~isempty(A.zbottom), ...
               'Provide alti, zmembrane, zbottom.');

        % Inject user/ENV values into properties (remain NaN if not set)
        obj.fs        = obj.localFromEnvOrArg('RBR_FS',        A.fs);
        obj.nfft      = obj.localFromEnvOrArg('RBR_NFFT',      A.nfft);
        obj.minFreq   = obj.localFromEnvOrArg('RBR_MINFREQ',   A.minFreq);
        obj.maxFreq   = obj.localFromEnvOrArg('RBR_MAXFREQ',   A.maxFreq);
        obj.igCutoff  = obj.localFromEnvOrArg('RBR_IGCUTOFF',  A.igCutoff);
        obj.windowSec = obj.localFromEnvOrArg('RBR_WINDOWSEC', A.windowSec);
        obj.critWL    = obj.localFromEnvOrArg('RBR_CRITWL',    A.critWL);

        % Require that all are now defined
        missing = {'fs','nfft','minFreq','maxFreq','igCutoff','windowSec','critWL'};
        missing = missing(cellfun(@(f) isnan(obj.(f)), missing));
        assert(isempty(missing), 'Missing required parameters: %s', strjoin(missing, ', '));

        % Geometry
        hd = A.zmembrane - A.zbottom; % sensor height above bed

        % ---- 1) RBR .rsk ----
        [t_rbr_num, p_raw_Pa] = obj.readRSKpressure(A.RSK);    % Pa

        % ---- 2) Meteo CSV ----
        [t_met_num, p_atm_kPa, T_C] = obj.readMeteoCSV(A.MeteoCSV);
        p_atm_Pa = obj.fillmissing_linear(p_atm_kPa*1000);     % kPa->Pa + fill gaps
        p_atm_interp = interp1(t_met_num, p_atm_Pa, t_rbr_num, 'spline', 'extrap');
        T_interp     = interp1(t_met_num, T_C,      t_rbr_num, 'spline', 'extrap');

        hs = (obj.R_gas*(T_interp + 273.15)) ./ (obj.M_air * obj.g); % scale height [m]
        p_atm_corr = p_atm_interp .* exp(A.alti ./ hs);               % adjusted to MSL

        % Gauge pressure
        p_gauge = p_raw_Pa - p_atm_corr;                              % Pa

        % ---- 3)  spectra & bulk parameters ----
        ndelay  = round(obj.windowSec * obj.fs);
        nBlocks = floor(numel(p_gauge) / ndelay);

        Hs = nan(nBlocks,1); Hs_IG = Hs; Hs_SW = Hs;
        Tp = Hs; Tm01 = Hs; Tm02 = Hs;
        TimeB = nan(nBlocks,1); ABS_WL = TimeB; WL_CG = TimeB;

        for ii = 1:nBlocks
            idx  = (1:ndelay) + (ii-1)*ndelay;
            segP = p_gauge(idx);
            tt   = t_rbr_num(idx);

            TimeB(ii) = mean(tt);
            pm        = mean(segP);
            ABS_WL(ii)= pm/(obj.rho*obj.g) + hd;      % absolute WL 
            WL_CG(ii) = ABS_WL(ii) + A.zbottom;

            if ABS_WL(ii) <= obj.critWL
                continue;  % instrument likely emerged
            end

            % Pressure PSD 
            [ff, df, PP, ~, ~] = rbr.signal.spectrum_fabrice(segP, obj.fs, obj.nfft, ABS_WL(ii));

            % Frequency band selection
            m = (ff >= obj.minFreq) & (ff <= obj.maxFreq);
            ff = ff(m); PP = PP(m);
            if isempty(ff), continue; end

            % Hydrostatic conversion: pressure -> elevation PSD
            PP = PP ./ (obj.rho*obj.g).^2;           % [m^2/Hz]

            % Bulk wave parameters
            [Hs(ii), Tp(ii), Tm01(ii), Tm02(ii), Hs_IG(ii), Hs_SW(ii)] = ...
                obj.bulkParams(ff, PP, df, obj.igCutoff);
        end

        % ---- pack outputs ----
        out = struct();
        out.Time            = TimeB;
        out.ABSOLUTE_WL     = ABS_WL;
        out.WL_CGVD2013     = WL_CG;
        out.spec.Hs         = Hs;
        out.spec.Hs_IG      = Hs_IG;
        out.spec.Hs_SW      = Hs_SW;
        out.spec.Tp         = Tp;
        out.spec.Tm01       = Tm01;
        out.spec.Tm02       = Tm02;
        out.meta.fs         = obj.fs;
        out.meta.nfft       = obj.nfft;
        out.meta.minFreq    = obj.minFreq;
        out.meta.maxFreq    = obj.maxFreq;
        out.meta.igCutoff   = obj.igCutoff;
        out.meta.windowSec  = obj.windowSec;
        out.meta.critWL     = obj.critWL;
        out.meta.hd         = hd;
        out.meta.alti       = A.alti;
        out.meta.zmembrane  = A.zmembrane;
        out.meta.zbottom    = A.zbottom;
    end
end

methods (Access = private)
    function [t_num, p_Pa] = readRSKpressure(~, rskFile)
        % Read RBR .rsk with RSKtools.
        RSK  = RSKopen(rskFile);
        Data = RSKreaddata(RSK);
        t    = Data.data.tstamp;
        if isdatetime(t), t_num = datenum(t); else, t_num = t; end
        v    = Data.data.values;
        p_Pa = v(:,1) * 1e4;      % dbar -> Pa
    end

    function [t_num, p_kPa, T_C] = readMeteoCSV(~, csvPath)
        A = readmatrix(csvPath);
        assert(size(A,2) >= 3, 'Meteo CSV must have 3 columns (Time, Press, Temp).');
        t_num = A(:,1); p_kPa = A(:,2); T_C = A(:,3);
    end

    function x = fillmissing_linear(~, x)
        idx = isnan(x);
        if any(idx)
            xi = find(~idx);
            x  = interp1(xi, x(~idx), 1:numel(x), 'linear', 'extrap').';
        end
    end

    function [Hs, Tp, Tm01, Tm02, HsIG, HsSW] = bulkParams(~, ff, PP, df, igCut)
        m0 = sum(PP)        * df;
        m1 = sum(ff.*PP)    * df;
        m2 = sum(ff.^2.*PP) * df;
        Hs   = 4*sqrt(m0);
        [~,pk] = max(PP);   Tp = 1/ff(pk);
        Tm01 = m0/max(m1,eps);
        Tm02 = sqrt(m0/max(m2,eps));
        Iig  = ff < igCut;
        HsIG = 4*sqrt(sum(PP(Iig))*df);
        HsSW = 4*sqrt(sum(PP(~Iig))*df);
    end

    function val = localFromEnvOrArg(~, envName, argVal)
        % If argVal is empty/NaN, try environment variable (string) then str2double
        if ~isempty(argVal), val = argVal; return; end
        s = getenv(envName);
        if isempty(s)
            val = NaN;
        else
            tmp = str2double(s);
            val = tmp;
            if isnan(val) && ~isempty(s)
                % allow non-numeric env for some future cases
                val = NaN;
            end
        end
    end
end
end



    
