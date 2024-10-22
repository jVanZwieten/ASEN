clear
close all
addpath('..')

registerA = ones(1, 4);
codeA = codeFromRegister(registerA, [4 3], 15)
codeB = codeFromRegister(registerA, [4 1], 15)

R_A = correlate(codeA, codeA)
R_B = correlate(codeA, codeB)

function code = codeFromRegister(register, keys, iterations)
    code = zeros(1, iterations);
    for i = 1:iterations
        code(i) = registerOutput(register(end));
        register = [xor(register(keys(1)), register(keys(2))) register(1:end - 1)];
    end
end

function output = registerOutput(registerInput)
    output = registerInput*-2 + 1;
end

function R = correlate(codeA, codeB)
    assert(length(codeA) == length(codeB))
    R = NaN(1, length(codeA));
    for n = 1:length(R)
        codeB_n = Utilities.leftShift(codeB, n);
        R(n) = sum(codeA.*codeB_n);
    end
end