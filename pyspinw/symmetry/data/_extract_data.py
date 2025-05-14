""" Convert the raw data scraped from the Bilbo server to a more useful form """

tokens = []

with open("spacegroup_data_raw_bilbao.txt", 'r') as file:

    for line in file:
        tokens += [part.strip() for part in line.split("\t") if len(part) > 0]

# It's well ordered, so we can just do this
groups = [token for i, token in enumerate(tokens) if i % 2 == 1 ]

with open("spacegroup_names.txt", 'w') as file:
    for group in groups:
        file.write(f"{group}\n")

