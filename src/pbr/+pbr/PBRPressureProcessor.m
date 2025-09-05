@Date    : 2025-09-02
@Author  : Faten Zouaghi (faten_zouaghi@uqar.ca)

classdef PBRPressureBurstProcessor
% PBRPressureBurstProcessor
% Wave processing for DIY PBR pressure sensors:


    properties (Constant, Access=private)
        rho   = 1023;    % seawater density [kg/m^3]
        g     = 9.81;    % gravity [m/s^2]
        Rgas  = 8.314;   % universal gas constant [J/(mol*K)]
        Mair  = 0.02896; % dry-air molar mass [kg/mol]
    end

    % ---- required numerics (no defaults: kept NaN until user/ENV provides) ----
    properties
        fs        (1,1) double = NaN; 
        nfft      (1,1) double = NaN;  
        minFreq   (1,1) double = NaN;  
        maxFreq   (1,1) double = NaN;  
        igCutoff  (1,1) double = NaN;  
        segSec    (1,1) double = NaN;  
        stepMin   (1,1) double = NaN;  
        critWL    (1,1) double = NaN;  
    end

    methods
        function out = run(obj, varargin)
            % RUN  Execute the full PBR pipeline 
            %
            % Required (name–value OR ENV):
            %   'PBRDir'     or env PBR_DIR        : directory with *.txt
            %   'MeteoCSV'   or env PBR_METEO      : 
            %   'alti'       or env PBR_ALTI       : station altitude above MSL [m]
            %   'zmembrane'  or env PBR_ZMEMBRANE  : sensor membrane Z [m]
            %   'zbottom'    or env PBR_ZBOTTOM    : bottom Z [m]
            %   'fs'         or env PBR_FS         : Hz
            %   'nfft'       or env PBR_NFFT
            %   'minFreq'    or env PBR_MINFREQ
            %   'maxFreq'    or env PBR_MAXFREQ
            %   'igCutoff'   or env PBR_IGCUTOFF
            %   'segSec'     or env PBR_SEGSEC
            %   'stepMin'    or env PBR_STEPMIN
            %   'critWL'     or env PBR_CRITWL
            %
            % Notes:
    

            % ---------------- parse args ----------------
            ip = inputParser; ip.FunctionName = 'PBRPressureBurstProcessor.run';
            addParameter(ip,'PBRDir','',@(s)ischar(s)||isstring(s));
            addParameter(ip,'MeteoMAT','',@(s)ischar(s)||isstring(s));
            addParameter(ip,'alti',[],@isnumeric);
            addParameter(ip,'zmembrane',[],@isnumeric);
            addParameter(ip,'zbottom',[],@isnumeric);
            addParameter(ip,'fs',[],@isnumeric);
            addParameter(ip,'nfft',[],@isnumeric);
            addParameter(ip,'minFreq',[],@isnumeric);
            addParameter(ip,'maxFreq',[],@isnumeric);
            addParameter(ip,'igCutoff',[],@isnumeric);
            addParameter(ip,'segSec',[],@isnumeric);
            addParameter(ip,'stepMin',[],@isnumeric);
            addParameter(ip,'critWL',[],@isnumeric);
            parse(ip,varargin{:});
            A = ip.Results;

            % ---------------- resolve inputs via ENV if missing ----------------
            pbrDir    = string( obj.fromArgOrEnv(A.PBRDir,   'PBR_DIR') );
            meteoPath = string( obj.fromArgOrEnv(A.MeteoMAT, 'PBR_METEO') );

            alti       = obj.numFromArgOrEnv(A.alti,      'PBR_ALTI');
            zmembrane  = obj.numFromArgOrEnv(A.zmembrane, 'PBR_ZMEMBRANE');
            zbottom    = obj.numFromArgOrEnv(A.zbottom,   'PBR_ZBOTTOM');

            obj.fs       = obj.numFromArgOrEnv(A.fs,       'PBR_FS');
            obj.nfft     = obj.numFromArgOrEnv(A.nfft,     'PBR_NFFT');
            obj.minFreq  = obj.numFromArgOrEnv(A.minFreq,  'PBR_MINFREQ');
            obj.maxFreq  = obj.numFromArgOrEnv(A.maxFreq,  'PBR_MAXFREQ');
            obj.igCutoff = obj.numFromArgOrEnv(A.igCutoff, 'PBR_IGCUTOFF');
            obj.segSec   = obj.numFromArgOrEnv(A.segSec,   'PBR_SEGSEC');
            obj.stepMin  = obj.numFromArgOrEnv(A.stepMin,  'PBR_STEPMIN');
            obj.critWL   = obj.numFromArgOrEnv(A.critWL,   'PBR_CRITWL');

            % ---------------- fail-fast checks ----------------
            assert(~isempty(pbrDir)    && isfolder(pbrDir),    'PBR directory not found.');
            assert(~isempty(meteoPath) && isfile(meteoPath),   'Meteo CSV not found.');

            obj.requireDefined({'fs','nfft','minFreq','maxFreq','igCutoff','segSec','stepMin','critWL'});
            assert(~isempty(alti) && ~isempty(zmembrane) && ~isempty(zbottom), ...
                   'Provide alti, zmembrane, zbottom.');

            % geometry
            hd = zmembrane - zbottom;   % sensor height above bed

            % ---------------- 1) Read PBR *.txt → time & pressure ----------------
            [t_sensor_num, pressure_Pa] = obj.readPBRdir(pbrDir);

            % ---------------- 2) Meteo MAT & barometric leveling ----------------
            [t_met_num, p_kPa, T_C] = obj.readMeteoCSV(meteoPath);
            p_Pa = obj.fillmissing_linear(p_kPa*1000);  % kPa→Pa + linear fill 
            p_atm   = interp1(t_met_num, p_Pa, t_sensor_num, 'spline', 'extrap');
            T_interp= interp1(t_met_num, T_C,  t_sensor_num, 'spline', 'extrap');

            % isothermal scale height
            hs = (obj.Rgas*(T_interp + 273.15)) ./ (obj.Mair*obj.g);
            p_atm_corr = p_atm .* exp(alti ./ hs);  

            % gauge pressure
            Pressure = pressure_Pa - p_atm_corr;

            % ---------------- 3) Burst loop ----------------
            ndelay    = round(obj.segSec * obj.fs);
            t0        = t_sensor_num(1);
            t1        = t_sensor_num(end);
            step_days = (obj.stepMin/60/24);
            starts    = t0 : step_days : t1;
            nB        = numel(starts);

            Hs = nan(nB,1); Tp = Hs; Tm01 = Hs; Tm02 = Hs; Hs_IG = Hs; Hs_SW = Hs;
            TimeB = nan(nB,1); ABS_WL = TimeB; WL_CG = TimeB;

            for i=1:nB
                [~, sidx] = min(abs(t_sensor_num - starts(i)));
                eidx = sidx + ndelay - 1;
                if eidx > numel(Pressure), break; end

                segP = Pressure(sidx:eidx);
                segT = t_sensor_num(sidx:eidx);
                TimeB(i) = mean(segT);

                pm = mean(segP);
                ABS_WL(i) = pm/(obj.rho*obj.g) + hd;
                WL_CG(i)  = ABS_WL(i) + zbottom;

                if ABS_WL(i) <= obj.critWL
                    continue; % likely emerged
                end

                % spectra
                [ff, df, PPp, ~, ~] = rbr.signal.spectrum_fabrice(segP, obj.fs, obj.nfft, ABS_WL(i));

                % select band
                K = (ff >= obj.minFreq) & (ff <= obj.maxFreq);
                ff = ff(K);  PPp = PPp(K);
                if isempty(ff), continue; end

                %  pressure → elevation PSD
                PP = PPp ./ (obj.rho*obj.g).^2;

                % bulk wave params (+ IG / sea–swell split)
                [Hs(i), Tp(i), Tm01(i), Tm02(i), Hs_IG(i), Hs_SW(i)] = ...
                    obj.bulkParams(ff, PP, df, obj.igCutoff);
            end

            % ---------------- 4) pack outputs ----------------
            out = struct();
            out.Time            = TimeB;             % datenum
            out.ABSOLUTE_WL     = ABS_WL;           % m
            out.WL_CGVD2013     = WL_CG;            % m
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
            out.meta.segSec     = obj.segSec;
            out.meta.stepMin    = obj.stepMin;
            out.meta.critWL     = obj.critWL;
            out.meta.hd         = hd;
            out.meta.alti       = alti;
            out.meta.zmembrane  = zmembrane;
            out.meta.zbottom    = zbottom;
        end
    end

    methods (Access=private)
        function [t_num, p_Pa] = readPBRdir(~, pbrDir)

            L = dir(fullfile(pbrDir, '*.txt'));
            assert(~isempty(L), 'No *.txt found in PBRDir.');

            data_DIY = [];
            for k = 1:numel(L)
                T = readtable(fullfile(pbrDir, L(k).name), 'Delimiter',' ', 'ReadVariableNames', false);
                data_DIY = [data_DIY; T]; %#ok<AGROW>
            end

            start_burst = find(data_DIY{:,1}==0);
            fmt = 'dd-MM-yyyyHH:mm:ss.SSS';

            % prealloc using total length
            t_sensor = NaT(height(data_DIY),1);  % timetable-compatible
            pressure = nan(height(data_DIY),1);

            for i = 1:numel(start_burst)-1
                i0 = start_burst(i);
                i1 = start_burst(i+1)-1;

                d  = num2str(data_DIY{i0,4});  % day
                m  = num2str(data_DIY{i0,5});  % month
                y  = num2str(data_DIY{i0,6});  % year
                hhmmss = char(data_DIY{i0,3}); % 'HH:mm:ss'

                t0 = datetime( strcat(d,'-',m,'-',y, hhmmss,'.',sprintf('%03d',0)) , ...
                               'InputFormat', fmt, 'TimeZone','UTC', 'Format', fmt);

                ms  = data_DIY{i0:i1,2};
                t_sensor(i0:i1) = t0 + milliseconds(ms);
                pressure(i0:i1) = data_DIY{i0:i1,8} * 100;  % hPa → Pa
            end

            % drop unfilled tails (last partial block, if any)
            K = ~isnat(t_sensor) & ~isnan(pressure);
            t_sensor = t_sensor(K);
            pressure = pressure(K);
        end

        function [t_num, p_kPa, T_C] = readMeteoMAT(~, metPath)
            % READMETEOMAT  Expect variables: Time_UTC (datenum or datetime), Press (kPa), Temperature (°C)
            S = load(metPath);
            % try to locate struct inside MAT (unnamed or named)
            fn = fieldnames(S);
            if numel(fn)==1 && isstruct(S.(fn{1})) && all(isfield(S.(fn{1}), {'Time_UTC','Press','Temperature'}))
                R = S.(fn{1});
            elseif all(isfield(S, {'Time_UTC','Press','Temperature'}))
                R = S;
            else
                error('Meteo MAT must contain Time_UTC, Press (kPa), Temperature (°C).');
            end

            Time_UTC = R.Time_UTC;
            if isdatetime(Time_UTC)
                t_num = datenum(Time_UTC);
            else
                t_num = Time_UTC;
            end
            p_kPa = R.Press(:);
            T_C   = R.Temperature(:);
        end

        function x = fillmissing_linear(~, x)
            % Linear gap-fill in index space (column vector expected)
            x = x(:);
            idx = isnan(x);
            if any(idx)
                xi  = find(~idx);
                x   = interp1(xi, x(~idx), (1:numel(x))', 'linear', 'extrap');
            end
        end

        function [Hs, Tp, Tm01, Tm02, HsIG, HsSW] = bulkParams(~, ff, PP, df, igCut)
            m0  = sum(PP)        * df;
            m1  = sum(ff.*PP)    * df;
            m2  = sum(ff.^2.*PP) * df;

            Hs   = 4*sqrt(m0);
            [~,pk] = max(PP);  Tp = 1/ff(pk);
            Tm01 = m0/max(m1,eps);
            Tm02 = sqrt(m0/max(m2,eps));

            Iig  = ff < igCut;
            HsIG = 4*sqrt(sum(PP(Iig))*df);
            HsSW = 4*sqrt(sum(PP(~Iig))*df);
        end

        function v = fromArgOrEnv(~, arg, envName)
            % Return char from arg or getenv; empty char if nothing.
            if ~(isstring(arg) || ischar(arg)) || strlength(string(arg))==0
                v = getenv(envName);
            else
                v = char(arg);
            end
        end

        function val = numFromArgOrEnv(~, arg, envName)
            % Numeric from arg or ENV; returns NaN if unresolved/non-numeric.
            if ~isempty(arg)
                val = arg;
                return;
            end
            s = getenv(envName);
            if isempty(s)
                val = NaN;
            else
                val = str2double(s);
                if isnan(val), val = NaN; end
            end
        end

        function requireDefined(obj, names)
            % Ensure all listed properties are numerically defined (not NaN).
            for k=1:numel(names)
                if isnan(obj.(names{k}))
                    error('Missing required numeric parameter: %s', names{k});
                end
            end
        end
    end
end
