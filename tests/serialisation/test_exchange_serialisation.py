import pytest
from pyspinw import LatticeSite

import numpy as np

from pyspinw.exchange import all_exchanges, Exchange

sites = [
    LatticeSite(0,0,1/2, 0,0,1),
    LatticeSite(1/2,0,1/2, 2,0,1),
    LatticeSite(0,1/2,1/2, 0,4,1),
    LatticeSite(0,1/2,1/2, 0,0,5),
    LatticeSite(1/2,1/2,1/2, 0,1,1),
]


site_pairs = [(site1, site2)
              for index, site1 in enumerate(sites)
              for site2 in sites[:index]]

rng = np.random.default_rng(42069)

@pytest.mark.parametrize("pair", site_pairs)
@pytest.mark.parametrize("exchange_class", all_exchanges)
def test_exchanges_serialise(pair: tuple[LatticeSite, LatticeSite], exchange_class: type[Exchange]):
    parameters = {parameter: rng.random() for parameter in exchange_class.parameters}

    exchange = exchange_class(*pair, **parameters)

    json = exchange.serialise()

    deserialised = Exchange.deserialise(json)

    assert isinstance(deserialised, exchange_class)

    assert np.all(exchange.site_1.ijk == deserialised.site_1.ijk)
    assert np.all(exchange.site_2.ijk == deserialised.site_2.ijk)

    assert np.all(exchange.site_1.spin_data == deserialised.site_1.spin_data)
    assert np.all(exchange.site_2.spin_data == deserialised.site_2.spin_data)

    print(parameters)
    print(deserialised.__dict__)

    for parameter in parameters:
        assert deserialised.__dict__["_" + parameter] == parameters[parameter]