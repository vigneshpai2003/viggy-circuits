class Switch:
    def __init__(self):
        self.__closed = True

    @property
    def isClosed(self) -> bool:
        return self.__closed

    @property
    def isOpen(self) -> bool:
        return not self.isClosed

    def open(self):
        self.__closed = False

    def close(self):
        self.__closed = True

    def toggle(self):
        self.__closed = not self.__closed
