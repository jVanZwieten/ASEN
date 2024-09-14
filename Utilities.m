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
    end
end