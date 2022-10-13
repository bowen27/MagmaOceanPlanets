# test_with_pytest.py
import numpy as np
from mop import parameters, calc_test
par = calc_test()

def test_parameters():
    from mop import parameters
    x = parameters()

    # Assert parameters are all positive
    assert x.rp > 0
    assert x.ndeg > 0
    assert x.ddeg > 0 

    assert x.dt > 0
    assert x.tmin >= 0
    assert x.tmax > 0
    assert x.tmin < x.tmax, "max time < min time"
    assert x.dt/x.dx < 0.4, "CFL Criterion not met"


def test_get_para():
    from mop import get_para
    # Didn't check sinusoidal case

    # test default forcing
    Fnet_array, p0_array = get_para(par, option=0)
    assert not np.any(Fnet_array), "Default forcing is non-zero"
    assert all(elem == Fnet_array[0] for elem in Fnet_array), 'Default forcing is non uniform'

    # test option = 1 forcing
    Fnet_array, p0_array = get_para(par, option=1)
    assert all(elem == Fnet_array[0] for elem in Fnet_array), 'Option = 1 (Uniform forcing) forcing is non uniform'

    # test option = random number that is not specified
    Fnet_array, p0_array = get_para(par, option=100)
    assert not np.any(Fnet_array), "Default forcing is non-zero"
    assert all(elem == Fnet_array[0] for elem in Fnet_array), 'Default forcing is non uniform'


def test_get_esat():
    from mop import get_esat
    pass

    # Maybe could create a dependent class of test values 

def test_get_Cm():
    pass

    


