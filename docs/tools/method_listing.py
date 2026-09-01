""" Makes files used for helping to write cross class documentation"""

import pyspinw
import inspect

file = open("classes_and_methods.txt", 'w')
padding = 50

def printy(*args):
    string = " ".join(str(x) for x in args)
    print(string)
    file.write(string)
    file.write(" "*(padding-len(string)))
    file.write("\n")


for name in pyspinw.__all__:
    obj = pyspinw.__dict__[name]

    if inspect.isfunction(obj) and not inspect.isclass(obj):
        printy("Function:", name)
        continue

    if inspect.isclass(obj):
        printy("Class:", name)

        methods = inspect.getmembers(obj, predicate=inspect.isfunction)
        properties = inspect.getmembers(obj, predicate=lambda x: isinstance(x, property))

        for method_name, method_obj in methods:
            if not method_name.startswith("_"):
                printy("  Method:", method_name)

        for property_name, property_obj in properties:
            printy("  Property:", property_name)


file.close()