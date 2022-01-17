from circtester import CircuitRunner

def test_circrunner_getsso():
    dist1={"1111": .9, "0000": 0.1}
    dist2={"1111": .9, "0000": 0.1}
    assert CircuitRunner.get_sso(dist1, dist2) == 1.0