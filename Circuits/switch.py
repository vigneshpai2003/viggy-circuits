class Switch:
    def __init__(self):
        self.__closed = True

    @property
    def is_closed(self) -> bool:
        return self.__closed

    @property
    def is_open(self) -> bool:
        return not self.is_closed

    def open(self):
        self.__closed = False

    def close(self):
        self.__closed = True

    def toggle(self):
        self.__closed = not self.__closed
