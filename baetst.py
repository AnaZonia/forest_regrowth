class Person:
    def __init__(self, name):
        self.name = name

    def greet(self):
        return f"Hello, my name is {self.name}."


class bf(Person):  # bf class inherits from Person class
    # init makes qualifier things for the class: every bf needs a gf
    def __init__(self, name, gf):
        super().__init__(name)  # Call the init method of the Person class instead of rewriting everything!
        self.gf = gf

    def love(self):
        if self.gf != "Ana" and self.name == "Hossein":
            return "OH NO WHERE IS MAH BAE."
        return f"Hello, my name is {self.name} and {self.gf} is MAH BAE."


def main():
    # Create two Person objects
    person1 = Person("Alice")
    person2 = Person("Bob")
    Hossein = Person("Hossein")
    Ana = Person("Ana")
    case1 = bf(Hossein.name, Ana.name)
    case2 = bf("Juan", "Saba")
    case3 = bf("Hossein", "Saba")

    # Call the greet method for both objects
    print(person1.greet())
    print(person2.greet())
    print(case1.love())
    print(case2.love())
    print(case3.love())


if __name__ == "__main__":
    main()