"""Tiny FreeSASA test double for subprocess integration tests."""

__version__ = "test-freesasa"
ShrakeRupley = "ShrakeRupley"


class Parameters:
    def __init__(self, values):
        self.values = dict(values)


class _Result:
    def __init__(self, area):
        self._area = float(area)

    def totalArea(self):
        return self._area


def calcCoord(coords, radii, parameters):
    # Integration tests exercise pipeline wiring, not FreeSASA numerics.
    # A deterministic per-atom surface surrogate is sufficient here.
    return _Result(sum(4.0 * 3.141592653589793 * (r + 1.4) ** 2 for r in radii))
