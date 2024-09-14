classdef GnssUtilities
    methods(Static)
        function [G, G1, G2] = generateCA(phaseSelector, bits)
            G1 = ones(1, 10);
            G2 = ones(1, 10);
            [G, G1, G2] = GnssUtilities.generateCAFromRegisters(phaseSelector, G1, G2, bits);
        end
        
        function [G, G1, G2] = generateCAFromRegisters(phaseSelector, G1, G2, bits)
            G = NaN(1, bits);
            
            for i = 1:bits
                G1i = G1(10);
                
                S = [G2(phaseSelector(1)) G2(phaseSelector(2))];
                G2i = mod(sum(S), 2);
                
                XGi = mod(sum([G1i G2i]), 2);
                G(i) = XGi;
                
                G1Feedback = [G1(3) G1(10)];
                G1 = [mod(sum(G1Feedback), 2) G1(1:9)];
                
                G2Feedback = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
                G2 = [mod(sum(G2Feedback), 2) G2(1:9)];
            end
        end
        
        function Rn = autoCorrelateFunction(CA)
            Rn = [1:length(CA); ones(1, length(CA))];
            Rn(2, :) = arrayfun(@(x) GnssUtilities.autoCorrelateN(CA, x), Rn(1, :));
        end
        
        function R = autoCorrelateN(CA, n)
            x = GnssUtilities.CaToX(CA);
            xShiftN = Utilities.leftShift(x, n);
            R = GnssUtilities.correlate(x, xShiftN);
        end
        
        function R = correlate(x1, x2)
            R = sum(x1.*x2)/length(x1);
        end

        function Rn = crossCorrelateFunction(CA1, CA2)
            x1 = GnssUtilities.CaToX(CA1);
            x2 = GnssUtilities.CaToX(CA2);            
            Rn =  GnssUtilities.crossCorrelateFunctionX(x1, x2);
        end

        function Rn = crossCorrelateFunctionX(x1, x2)
            assert(length(x1)==length(x2))
            Rn = [1:length(x1); ones(1, length(x1))];
            Rn(2, :) = arrayfun(@(x) GnssUtilities.crossCorrelateNX(x1, x2, x), Rn(1, :));
        end

        function R = crossCorrelateN(CA1, CA2, n)
            x1 =  GnssUtilities.CaToX(CA1);
            x2 =  GnssUtilities.CaToX(CA2);
            R = GnssUtilities.crossCorrelateNX(x1, x2, n);
        end

        function R = crossCorrelateNX(x1, x2, n)
            x2 = Utilities.leftShift(x2, n);
            R = GnssUtilities.correlate(x1, x2);
        end

        function x = CaToX(CA)
            x = (CA==0)*1 + (CA==1)*-1;
        end
    end
end