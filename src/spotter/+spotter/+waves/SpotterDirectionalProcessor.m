classdef SpotterDirectionalProcessor
% Directional wave parameters from Spotter displacement CSV.
%
    properties (Access=private, Constant)
        DEF = struct('Fs',2.5,'SegLength',3600,'Bandpass',[0.05 0.5], ...
                     'Ncvec',8,'CI',95,'Window',"hanning",'P',0);
    end

    methods
        function out = run(obj, varargin)
            % Parse des paramètres (name-value) + fallback env vars
            p = inputParser;   %#ok<*TIPD> 
            p.FunctionName = 'SpotterDirectionalProcessor.run';
            addParameter(p,'CSV',   "",           @(s)isstring(s) || ischar(s));
            addParameter(p,'Fs',    obj.DEF.Fs,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
            addParameter(p,'SegLength', obj.DEF.SegLength, @(x)isnumeric(x)&&isscalar(x)&&x>0);
            addParameter(p,'Bandpass', obj.DEF.Bandpass, @(v)isnumeric(v)&&numel(v)==2&&all(v>0));
            addParameter(p,'Ncvec', obj.DEF.Ncvec, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
            addParameter(p,'CI',    obj.DEF.CI,   @(x)isnumeric(x)&&isscalar(x)&&x>0 && x<=100);
            addParameter(p,'Window',obj.DEF.Window,@(s)isstring(s)||ischar(s));
            addParameter(p,'P',     obj.DEF.P,    @(x)isnumeric(x)&&isscalar(x));
            parse(p, varargin{:});
            args = p.Results;

            % CSV local : arg > env var > erreur
            csvPath = string(args.CSV);
            if (csvPath=="") || ~isfile(csvPath)
                envCSV = string(getenv('SPOTTER_CSV'));
                if envCSV~="" && isfile(envCSV)
                    csvPath = envCSV;
                else
                    error('CSV file not provided or not found. Pass ''CSV'',<path> or setenv(''SPOTTER_CSV'',<path>).');
                end
            end

            % Dépendances requises
            assert(exist('psd_calc','file')==2, 'Missing dependency: psd_calc.m not on MATLAB path.');
            assert(exist('csd_calc','file')==2, 'Missing dependency: csd_calc.m not on MATLAB path.');

            % 1) Lecture CSV + régularisation 2.5 Hz
            [temps, disp_est, disp_nord, disp_vertical] = obj.readAndRegularize(csvPath, args.Fs);

            % 2) Filtre passe-bande (Butterworth) + filtfilt (zéro phase)
            [disp_est, disp_nord, disp_vertical] = obj.bandpassAll(disp_est, disp_nord, disp_vertical, args.Fs, args.Bandpass);

            % 3) Segmentation -> calcul PSD/CSD par segment
            [tm, ff, Szz, Sxx, Syy, Sxy, Qyz, Qxz] = obj.blockSpectra(temps, disp_est, disp_nord, disp_vertical, ...
                args.Fs, args.SegLength, args.Ncvec, args.CI, string(args.Window), args.P);

            % 4) Sélection bande gravitaire
            J = (ff>0.05) & (ff<0.5);
            ff  = ff(J);
            Szz = Szz(:,J); Sxx = Sxx(:,J); Syy = Syy(:,J);
            Sxy = Sxy(:,J); Qyz = Qyz(:,J); Qxz = Qxz(:,J);

            % 5) Coeffs directionnels
            den = sqrt( Szz .* (Sxx + Syy) );
            a1 = Qxz ./ den;
            b1 = Qyz ./ den;
            a2 = (Sxx - Syy) ./ (Sxx + Syy);
            b2 = (2 .* Sxy) ./ (Sxx + Syy);

            % 6) Moments + paramètres bulk
            dff = gradient(ff);
            m0 = obj.rowSum(Szz .* dff);
            m1 = obj.rowSum(Szz .* (ff.*dff));
            m2 = obj.rowSum(Szz .* ((ff.^2).*dff));

            Hs  = 4.*sqrt(m0);
            [~,pk] = max(Szz,[],2);
            fp  = ff(pk);   Tp = 1./fp(:);
            Tm01 = m0 ./ m1;
            Tm02 = sqrt(m0 ./ m2);

            % Moyennes spectrales (Kuik-style)
            a1mean = obj.spectralMean(Szz, a1, dff, m0);
            b1mean = obj.spectralMean(Szz, b1, dff, m0);
            a2mean = obj.spectralMean(Szz, a2, dff, m0);
            b2mean = obj.spectralMean(Szz, b2, dff, m0);

            % Directions et spreads (renommages demandés)
            Dm = mod( 270 - (180/pi).*atan2(b1mean, a1mean), 360 );
            ai = a1(sub2ind(size(a1),(1:numel(pk))',pk));
            bi = b1(sub2ind(size(b1),(1:numel(pk))',pk));
            Dp = mod( 270 - (180/pi).*atan2(bi, ai), 360 );

            MDS = mod( (180/pi).*sqrt( 2.*(1 - sqrt(a1mean.^2 + b1mean.^2)) ), 360 );
            PDS = mod( (180/pi).*sqrt( 2.*(1 - sqrt(ai.^2 + bi.^2)) ), 360 );

            % Sortie (structure)
            out = struct();
            out.tm   = tm(:);
            out.ff   = ff(:)';
            out.Szz  = Szz; out.Sxx=Sxx; out.Syy=Syy; out.Sxy=Sxy; out.Qyz=Qyz; out.Qxz=Qxz;
            out.Hs   = Hs(:); out.Tp=Tp(:); out.Tm01=Tm01(:); out.Tm02=Tm02(:);
            out.a1mean=a1mean(:); out.b1mean=b1mean(:); out.a2mean=a2mean(:); out.b2mean=b2mean(:);
            out.Dm   = Dm(:); out.Dp=Dp(:); out.MDS=MDS(:); out.PDS=PDS(:);
        end
    end

    methods (Access=private)
        function [temps_num, de, dn, dz] = readAndRegularize(~, csvPath, Fs)
            T = readtable(csvPath);
            A = table2array(T);
            yr=A(:,1); mo=A(:,2); da=A(:,3); hh=A(:,4); mm=A(:,5); ss=A(:,6); ms=A(:,7);
            tt = datetime(yr,mo,da,hh,mm,ss+ms/1000,'TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss.SSS');

            de  = A(:,8);   % est (X)
            dn  = A(:,9);   % nord (Y)
            dz  = A(:,10);  % vertical (Z)

            M = timetable(tt,dz,dn,de);
            M = retime(M,'regular','linear','TimeStep',seconds(1/Fs));
            temps_num = datenum(M.tt);
            dz = M.dz; dn = M.dn; de = M.de;
        end

        function [de,dn,dz] = bandpassAll(~, de, dn, dz, Fs, bp)
            fN = Fs/2;
            [B,A] = butter(2, bp./fN, 'bandpass');
            dz = filtfilt(B,A,dz);
            de = filtfilt(B,A,de);
            dn = filtfilt(B,A,dn);
        end

        function [tm, ff, Szz, Sxx, Syy, Sxy, Qyz, Qxz] = blockSpectra(~, temps, de, dn, dz, Fs, segsec, nc, ci, win, p)
            ndelay    = round(segsec * Fs);
            nbspectre = floor(numel(dz)/ndelay);

            Szz = []; Sxx = []; Syy = [];
            Sxy = []; Qyz = []; Qxz = [];
            ff  = []; tm  = zeros(nbspectre,1);

            for i=1:nbspectre
                idx = (1:ndelay) + (i-1)*ndelay;
                vz  = dz(idx);
                ve  = de(idx);
                vn  = dn(idx);
                tm(i) = mean(temps(idx));

                [psd_z, fHz, ~, ~] = spotter.signal.psd_calc(vz, 1/Fs, nc, ci, win, p);
                [psd_x, ~,   ~, ~] =spotter.signal.psd_calc(ve, 1/Fs, nc, ci, win, p);
                [psd_y, ~,   ~, ~] = spotter.signal.psd_calc(vn, 1/Fs, nc, ci, win, p);

                [csd_xy, ~, ~] = spotter.signal.csd_calc(ve, vn, 1/Fs, nc, win, p);
                [csd_yz, ~, ~] =spotter.signal.csd_calc(vn, vz, 1/Fs, nc, win, p);
                [csd_xz, ~, ~] =spotter.signal.csd_calc(ve, vz, 1/Fs, nc, win, p);

                if isempty(ff), ff = fHz(:)'; end
                Szz(i,:) = psd_z(:)';  Sxx(i,:) = psd_x(:)';  Syy(i,:) = psd_y(:)';
                Sxy(i,:) = real(csd_xy(:)');    
                Qyz(i,:) = -imag(csd_yz(:)');  
                Qxz(i,:) = -imag(csd_xz(:)');  
            end
        end

        function s = rowSum(~, M)
            % Somme ligne par ligne d'une matrice 2D
            s = sum(M, 2);
        end

        function m = spectralMean(~, Szz, CC, dff, m0)
            % moyenne spectrale pondérée par Szz/m0
            n = size(Szz,1);
            m = zeros(n,1);
            for j=1:n
                m(j) = (1./m0(j)) * sum( Szz(j,:) .* CC(j,:) .* dff );
            end
        end
    end
end
