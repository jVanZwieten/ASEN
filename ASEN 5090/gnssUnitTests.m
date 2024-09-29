addpath('..')

range = 1500;
inputUserEcef = Utilities.radius_earth*[1; 0; 0];
inputSatelliteEcef = (Utilities.radius_earth + range)*[1; 0; 0];
expected = [0 pi/2 range];
computeAzElRangeTest(inputUserEcef, inputSatelliteEcef, expected)

function computeAzElRangeTest(userECEF, satECEF, expected)
    [az, el, range] = computeAzElRange(userECEF, satECEF);
    assert(az == expected(1))
    assert(el == expected(2))
    assert(range == expected(3))
end