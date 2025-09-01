classdef RBRPressureProcessor
% RBRPressureProcessor
% End-to-end pipeline: RBR pressure (.rsk) + meteo CSV -> wave time series
%
% Required external tool: RSKtools (RBR).
%   RSK = RSKopen('file.rsk'); RSK = RSKreaddata(RSK); etc.
%   (Install and add RSKtools to the MATLAB path before running.) 
%
% Outputs (per spectral block):
%   Time               : block-center time (datenum)
%   ABSOLUTE_WL        : absolute water level [m] (p/(rho*g) + hd)
%   WL_CGVD2013        : absolute WL relative to chart datum (ABSOLUTE_WL + zbottom)
%   Hs, Hs_IG, Hs_SW   : significant wave heights (total, infragravity, sea-swell)
%   Tp, Tm01, Tm02     : peak and mean wave periods
%
% Notes:
% - Meteo CSV must have 3 numeric columns: Time_UTC (MATLAB datenum), 
%   Press (kPa), Temperature (degC). No conversion to datetime is performed.
% - No dynamic pressure attenuation correction is applied; elevation PSD is
%   obtained by dividing pressure PSD by (rho*g)^2 (hydrostatic). 
%   You can later swap in a transfer function if desired.
%
% Author: (you)

properties (Constant)
    rho = 1023;                 % seawater density [kg/m^3]
    g   = 9.81;                 % gravity [m/s^2]
    R_gas = 8.314;              % [J/(mol*K)]
    M_air = 0.02896;            % [kg/mol]
end

properties
    fs = 4;                     % sampling frequency [Hz]
    nfft = 1024;                % FFT length
    minFreq = 0.0083;           % Hz
    maxFreq = 0.5;              % Hz
    igCutoff = 0.05;            % Hz (IG vs sea-swell split)
    windowSec = 20*60;          % one spectrum every 20 min
    critWL = 0.35;              % minimum water depth to accept a spectrum [m]
end

methods
    function out = run(obj, varargin)
        % RUN  Compute wave parameters from an RBR .rsk and a meteo CSV.
        %
        % Usage:
        %   proc = rbr.waves.RBRPressureProcessor;
        %   out  = proc.run('RSK','C:\data\file.rsk', ...
        %                   'MeteoCSV','C:\data\METEO.csv', ...
        %                   'alti',44.5,'zmembrane',-1.236,'zbottom',-1.52);
        %
        % Name-Value:
        %   'RSK'       : path to RBR .rsk
        %   'MeteoCSV'  : path to meteo CSV (Time_UTC, Press[kPa], Temperature[°C])
        %   'alti'      : met station altitude above MSL [m]
        %   'zmembrane' : sensor membrane elevation [m]
        %   'zbottom'   : bottom elevation [m]
        %   (optional overrides) 'fs','nfft','minFreq','maxFreq','igCutoff','windowSec','critWL'

        p = inputParser;
        addParameter(p,'RSK','',@(s)ischar(s)||isstring(s));
        addParameter(p,'MeteoCSV','',@(s)ischar(s)||isstring(s));
        addParameter(p,'alti',[],@isnumeric);
        addParameter(p,'zmembrane',[],@isnumeric);
        addParameter(p,'zbottom',[],@isnumeric);
        addParameter(p,'fs',obj.fs,@isnumeric);
        addParameter(p,'nfft',obj.nfft,@isnumeric);
        addParameter(p,'minFreq',obj.minFreq,@isnumeric);
        addParameter(p,'maxFreq',obj.maxFreq,@isnumeric);
        addParameter(p,'igCutoff',obj.igCutoff,@isnumeric);
        addParameter(p,'windowSec',obj.windowSec,@isnumeric);
        addParameter(p,'critWL',obj.critWL,@isnumeric);
        parse(p,varargin{:});
        args = p.Results;

        % allow env vars as a fallback (adds a tiny barrier for non-experts)
        if isempty(args.RSK),      args.RSK      = getenv('RBR_RSK');     end
        if isempty(args.MeteoCSV), args.MeteoCSV = getenv('RBR_METEO');   end
        assert(~isempty(args.RSK) && exist(args.RSK,'file')==2, 'RSK file not found.');
        assert(~isempty(args.MeteoCSV) && exist(args.MeteoCSV,'file')==2, 'Meteo CSV not found.');
        assert(~isempty(args.alti) && ~isempty(args.zmembrane) && ~isempty(args.zbottom), ...
            'Provide alti, zmembrane, zbottom.');

        % allow runtime overrides
        obj.fs       = args.fs;
        obj.nfft     = args.nfft;
        obj.minFreq  = args.minFreq;
        obj.maxFreq  = args.maxFreq;
        obj.igCutoff = args.igCutoff;
        obj.windowSec= args.windowSec;
        obj.critWL   = args.critWL;

        hd = args.zmembrane - args.zbottom;     % sensor height above bed

        % ---- 1) read RBR .rsk (absolute pressure) ----
        [t_rbr_num, p_raw_Pa] = obj.readRSKpressure(args.RSK);   % numeric datenum, Pa

        % ---- 2) read meteo CSV and barometric correction ----
        [t_met_num, p_atm_kPa, T_C] = obj.readMeteoCSV(args.MeteoCSV);
        p_atm_Pa = p_atm_kPa * 1000;                         % kPa -> Pa
        p_atm_Pa = obj.fillmissing_linear(p_atm_Pa);         % simple interp for gaps

        % interpolate meteo to RBR times (both are numeric datenums)
        p_atm_interp = interp1(t_met_num, p_atm_Pa, t_rbr_num, 'spline', 'extrap');
        T_interp     = interp1(t_met_num, T_C,      t_rbr_num, 'spline', 'extrap');

        % hypsometric adjustment to sea level (isothermal)
        hs = (obj.R_gas*(T_interp + 273.15)) ./ (obj.M_air * obj.g);    % [m]
        p_atm_corr = p_atm_interp .* exp(args.alti ./ hs);               % Pa

        % remove atmospheric pressure -> (quasi) gauge pressure
        p_gauge = p_raw_Pa - p_atm_corr;                                 % Pa

        % ---- 3) block processing & spectra ----
        ndelay    = round(obj.windowSec * obj.fs);
        nBlocks   = floor(numel(p_gauge) / ndelay);

        Hs    = nan(nBlocks,1);
        Tp    = nan(nBlocks,1);
        Tm01  = nan(nBlocks,1);
        Tm02  = nan(nBlocks,1);
        Hs_IG = nan(nBlocks,1);
        Hs_SW = nan(nBlocks,1);
        TimeB = nan(nBlocks,1);              % block-center time (datenum)
        ABS_WL= nan(nBlocks,1);              % absolute WL per block
        WL_CG = nan(nBlocks,1);

        for ii = 1:nBlocks
            idx  = (1:ndelay) + (ii-1)*ndelay;
            segP = p_gauge(idx);
            tt   = t_rbr_num(idx);
            TimeB(ii) = mean(tt);

            % mean absolute WL for the block (p/(rho*g) + hd)
            pm = mean(segP);
            ABS_WL(ii) = pm/(obj.rho*obj.g) + hd;
            WL_CG(ii)  = ABS_WL(ii) + args.zbottom;

            % skip if instrument emerges (threshold)
            if ABS_WL(ii) <= obj.critWL,  continue;  end

            % spectrum of pressure (Pa^2/Hz) using your routine
            % spectrum_fabrice MUST be available as rbr.signal.spectrum_fabrice
            [ff, df, PP, ~, ~] = rbr.signal.spectrum_fabrice(segP, obj.fs, obj.nfft, ABS_WL(ii));

            % band selection
            mask = (ff >= obj.minFreq) & (ff <= obj.maxFreq);
            ff = ff(mask);   PP = PP(mask);

            % pressure -> elevation PSD (hydrostatic)
            % (No dynamic transfer function; strictly divide by (rho*g)^2)
            PP = PP ./ (obj.rho*obj.g).^2;      % [m^2/Hz]

            if isempty(ff), continue; end

            % bulk parameters
            [Hs(ii), Tp(ii), Tm01(ii), Tm02(ii), Hs_IG(ii), Hs_SW(ii)] = ...
                obj.bulkParams(ff, PP, df, obj.igCutoff);
        end

        % ---- pack outputs ----
        out = struct();
        out.Time            = TimeB;                 % datenum
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
        out.meta.alti       = args.alti;
        out.meta.zmembrane  = args.zmembrane;
        out.meta.zbottom    = args.zbottom;
    end
end

methods (Access = private)

    function [t_num, p_Pa] = readRSKpressure(~, rskFile)
        % Read pressure time series from an RBR .rsk using RSKtools.
        RSK  = RSKopen(rskFile);
        Data = RSKreaddata(RSK);
        t    = Data.data.tstamp;
        if isdatetime(t), t_num = datenum(t); else, t_num = t; end
        v    = Data.data.values;
        p_Pa = v(:,1) * 1e4;         % RBR pressure channel in dbar -> Pa
    end

    function [t_num, p_kPa, T_C] = readMeteoCSV(~, csvPath)
        % Read meteo CSV (numeric matrix): [Time_UTC, Press(kPa), Temperature(C)]
        A = readmatrix(csvPath);
        assert(size(A,2)>=3, 'Meteo CSV must have 3 columns (Time, Press, Temp).');
        t_num = A(:,1);                 % already MATLAB serial time
        p_kPa = A(:,2);                 % kPa
        T_C   = A(:,3);                 % degC
    end

    function x = fillmissing_linear(~, x)
        % Fill NaNs by linear interpolation (endpoints held)
        idx = isnan(x);
        if any(idx)
            xi  = find(~idx);
            x   = interp1(xi, x(~idx), 1:numel(x), 'linear', 'extrap').';
        end
    end

    function [Hs, Tp, Tm01, Tm02, HsIG, HsSW] = bulkParams(~, ff, PP, df, igCut)
        % Spectral moments and bulk wave parameters from elevation PSD.
        m0 = sum(PP)       * df;
        m1 = sum(ff.*PP)   * df;
        m2 = sum(ff.^2.*PP)* df;

        Hs    = 4*sqrt(m0);
        [~,pk]= max(PP);            Tp = 1/ff(pk);
        Tm01  = m0/max(m1,eps);
        Tm02  = sqrt(m0/max(m2,eps));

        Iig   = ff < igCut;
        HsIG  = 4*sqrt(sum(PP(Iig))*df);
        HsSW  = 4*sqrt(sum(PP(~Iig))*df);
    end
end
end
