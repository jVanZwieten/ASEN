classdef Utilities
    properties(Constant)
        radius_earth = 6378137; % m, WGS-84 fundamental parameter
        flattening_earth = 1/298.257223563; % WGS-84 fundamental parameter
        eccentricity_earth = sqrt(2*Utilities.flattening_earth - Utilities.flattening_earth^2); % WGS-84 fundamental parameter
    end

    methods(Static)
        function hex = binVector2Hex(V)
            binary_str = num2str(V);
            binary_str = binary_str(~isspace(binary_str));
            hex = dec2hex(bin2dec(binary_str));
        end

        function clampedResult = clamp(x, upperBound, lowerBound)
            clampedResult = max(min(x, upperBound), lowerBound);
        end

        function s = days2s(d)
            s = d*24*60*60;
        end

        function varargout = decompose(V)
            for i = 1:size(V, 1)
                varargout{i} = V(i, :);
            end
        end

        function s = h2s(h)
            s = h*60*60;
        end

        function A = rightShift(A, n)
            A = A([(end - n + 1:end) (1:end - n)]);
        end

        function A = leftShift(A, n)
            A = A([(n + 1:end) (1:n)]);
        end

        function Ahat = UnitVector(A)
            Ahat = A/norm(A);
        end

        function multiplot(X, seriesLabels, axisLabels, figureTitle)
            n = size(X, 1) - 1;
            assert(length(seriesLabels) == n)
            assert(length(axisLabels) == n + 1)
            
            figure
            t = tiledlayout(n, 1);
            title(t, figureTitle)
            for i=2:(n + 1)
                nexttile
                plot(X(1, :), X(i, :))
                title(seriesLabels(i - 1))
                xlabel(axisLabels(1))
                ylabel(axisLabels(i))
            end
        end

        function multiplotY(Y, seriesLabels, axisLabels, figureTitle)
            n = size(Y, 1) - 1;
            assert(length(seriesLabels) == n)
            assert(length(axisLabels) == n + 1)
            
            figure
            t = tiledlayout(n, 1);
            title(t, figureTitle)
            for i=2:(n + 1)
                nexttile
                plot(Y(i, :), Y(1, :))
                title(seriesLabels(i - 1))
                xlabel(axisLabels(i))
                ylabel(axisLabels(1))
            end
        end

        function multiSeries(X, seriesLabels, axisLabels, figureTitle)
            figure
            hold on
            for i=2:size(X, 1)
                plot(X(1, :), X(i, :))
            end
            xlabel(axisLabels(1));
            ylabel(axisLabels(2));
            title(figureTitle);
            legend(seriesLabels)
        end

        function semicircles = rad2semicircle(rads)
            semicircles = rads/pi;
        end

        function d = s2days(s)
            d = s/60/60/24;
        end

        function rads = semicircle2rad(semicircles)
            rads = pi*semicircles;
        end

        function Norms = VectorizedNorms(varargin)
            Norms = zeros(size(varargin{1}));
            
            for i = 1:nargin
                Norms = Norms + varargin{i}.^2;
            end
            
            Norms = sqrt(Norms);
        end
    end
end