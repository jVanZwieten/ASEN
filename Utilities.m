classdef Utilities
    methods(Static)
        function hex = binVector2Hex(V)
            binary_str = num2str(V);
            binary_str = binary_str(~isspace(binary_str));
            hex = dec2hex(bin2dec(binary_str));
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
    end
end