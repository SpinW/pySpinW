""" Load in the tutorial files """

filenames = []
text = []
with open("tutorial_list.txt", 'r') as file:
    for line in file:
        parts = line.split(" ", maxsplit=1)
        filenames.append(parts[0]+".py")
        if len(parts) > 1:
            text.append(parts[1])
        else:
            text.append([])
