from .mass import Mass


class Ratio:
    def __init__(self, value):
        self._value = value

    @property
    def value(self):
        return self._value

    def __gt__(self, other):
        return self.value > other

    def __lt__(self, other):
        return self.value < other

    def __ge__(self, other):
        return self.value >= other

    def __le__(self, other):
        return self.value <= other

    def __eq__(self, other):
        return self.value == other


class HumidityRatio(Ratio):
    def __init__(self, solvant: Mass = None, solver: Mass = None, value: float = None):
        if value is None:
            value = solvant.kg / solver.kg

        super().__init__(value)
