def test_energy_exact():
    a_exact = 8 / (np.pi * 9)
    E_exact = -4 / (3 * np.pi)
    E_calc = energy([a_exact])
    assert abs(E_calc - E_exact) < 1e-6
