""" Load in the tutorial files """
import os.path

filenames = []
text = []

print(os.path.dirname(__file__))

path = os.path.abspath("../examples/tutorials/tutorial_list.txt")

with open(path, 'r') as file:
    for line in file:
        parts = line.split(" ", maxsplit=1)
        filenames.append(parts[0]+".py")
        if len(parts) > 1:
            text.append(parts[1])
        else:
            text.append([])
