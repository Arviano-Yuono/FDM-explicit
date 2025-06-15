class SingleNode:
    def __init__(self, name: str, x: float, T: float, type: str = "interior"):
        self.name = name
        self.x = x
        self.T = T
        self.type = type

    def __str__(self):
        return f"Node {self.name} at (x = {self.x}) with T={self.T} and type={self.type}"

    def __repr__(self):
        return self.__str__()