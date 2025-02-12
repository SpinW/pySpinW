from pyspinw.util.closures import closure
from pyspinw.util.group_generators import spglib_generators


def check_closed(number):
    print(f"{number}: ", end="")
    gens = spglib_generators(number)
    closed = closure(gens)

    if len(gens) != len(closed):
        print("Not closed")

        print("  Generators:")
        for gen in gens:
            print("   *", gen)

        print("  Closure:")
        for gen in closed:
            print("   *", gen)

        return False

    else:
        gens.sort()
        closed.sort()

        for a, b in zip(gens, closed):
            if a != b:
                print("Different")
                return False

        else:
            print("Closed")
            return True


all_closed = True
for i in range(1651):
    all_closed = check_closed(i+1) and all_closed

if all_closed:
    print("All closed")
else:
    print("Some groups were not closed")