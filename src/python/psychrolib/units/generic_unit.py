class GenericUnit:

    def __init__(self, base_value):
        self._base_value = base_value

    def __repr__(self):
        return f"{self._base_value} {self.__class__.__name__}"