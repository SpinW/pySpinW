from convert_tutorial_files import run as convert, tutorial_output_dir
from tutorial_list import filenames, text

with open(tutorial_output_dir / "index.md", 'w') as file:
    file.write("# Tutorials\n\n")

    for index, (filename, text) in enumerate(zip(filenames, text)):

        file.write(f" - [Tutorial {index+1}]({filename}/tutorial) {text}")

convert(filenames)