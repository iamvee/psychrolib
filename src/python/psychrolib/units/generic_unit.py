class GenericUnit:

    def __init__(self, base_value):
        self._base_value = base_value

    def __gt__(self, other):
        return self._base_value > other if isinstance(other, (float, int)) else other._base_value

    def __lt__(self, other):
        return self._base_value < other if isinstance(other, (float, int)) else other._base_value

    def __ge__(self, other):
        return self._base_value >= other if isinstance(other, (float, int)) else other._base_value

    def __le__(self, other):
        return self._base_value <= other if isinstance(other, (float, int)) else other._base_value

    
    def __eq__(self, other):
        return self._base_value == other if isinstance(other, (float, int)) else other._base_value

    def __repr__(self):
        return f"{self._base_value} {self.__class__.__name__}"