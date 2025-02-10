import numpy as np
from pyspinw.util.group_generators import Generator

def closure(generators: list[Generator], max_iters=1000) -> list[Generator]:
    """ Closure of magnetic space groups"""

    output_generators = [Generator(np.eye(3), np.zeros(3), 1, name="e")]

    for iter in range(max_iters):
        # print(f"Iteration {iter}")

        last_generators = output_generators.copy()
        for next_generator in generators:
            # apply this generator
            output_generators += [next_generator.and_then(generator) for generator in output_generators]

            # Put in unique order (uses Generator.__lt__)
            output_generators.sort()

            # Remove duplicates (uses Generator.__eq__)
            a = output_generators[0]
            new_output_generators = [a]
            for b in output_generators[1:]:
                if a != b:
                    new_output_generators.append(b)
                a = b

            output_generators = new_output_generators

        output_generators.sort()

        if len(last_generators) == len(output_generators):
            for a, b in zip(last_generators, output_generators):
                if a != b:
                    continue
            return output_generators

    else:
        raise Exception("Maximum iterations reached")

