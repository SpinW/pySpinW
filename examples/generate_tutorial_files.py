""" Generate files for tutorials that will go on the website """

import os
import shutil
from collections import defaultdict
from pathlib import Path
import runpy

input_dir = Path("tutorials")
output_dir = Path("tutorial_files")

spinw_tutorial_text = "\n_This tutorial mirrors MATLAB spinW Tutorial %s_\n\n"
md_placeholder = "%%%%% stdout "

# Clean directory
if os.path.exists(output_dir):

    for path in output_dir.iterdir():
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

else:
    os.mkdir(output_dir)

# Create files

for md_file in os.listdir(input_dir):
    if md_file.startswith("_"):
        continue

    if not md_file.endswith(".py"):
        continue

    base_name = md_file.split(".")[0]

    print(base_name)

    output_target_dir = output_dir / base_name

    os.mkdir(output_target_dir)

    # Parse file into blocks based on comment indentation

    comment_level = 0
    blocks = []
    current_block = []
    with open(input_dir / md_file, 'r') as input_file:
        for line in input_file:
            if line.startswith("## "):
                new_comment_level = 2
                add_line = line[3:]
            elif line.startswith("#"):
                new_comment_level = 1
                add_line = line[1:]
            else:
                new_comment_level = 0
                add_line = line


            if new_comment_level != comment_level:
                blocks.append((comment_level, current_block))
                current_block = []

            comment_level = new_comment_level
            current_block.append(add_line)

        blocks.append((comment_level, current_block))

    blocks = [(level, block) for level, block in blocks if not (level == 2 and len(block) == 0)]

    # for level, block in blocks:
    #     print(f"Level {level}:")
    #     for line in block:
    #         print("    ", line.rstrip())

    # Create the corresponding output files
    started = False
    skip_lines = 0
    stdout_capture_index = 0
    with open(output_target_dir / "tutorial.md", 'w') as md_file:
        with open(output_target_dir / "artifacts.py", 'w') as artifacts_file:
            artifacts_file.write("import sys\n"
                                 "import io\n"
                                 "old_stdout = sys.stdout\n"
                                 "buffer = io.StringIO()\n"
                                 "sys.stdout = buffer\n")

            for level, block in blocks:
                match level:

                    case 0: # It's code
                        if started:
                            # Skip any skipped lines
                            lines = block[skip_lines:]

                            # Clear any blank lines
                            lines = [line for line in lines if line.strip() != ""]

                            # Write to files
                            if len(lines) > 0:
                                md_file.write("```python\n")

                                for line in lines:
                                    artifacts_file.write(line)
                                    md_file.write(line)

                                md_file.write("```\n")

                        else:
                            pass #print(f"'{block}' comes before title tag, skipping")

                        skip_lines = 0

                    case 1: # It's text
                        if started:
                            for line in block[skip_lines:]:
                                md_file.write(line)
                        else:
                            pass # print(f"'{block}' comes before title tag, skipping")
                        skip_lines = 0

                    case 2: # It's a control tag
                        rest_is_code = False # Flag to treat the rest of this block as code to be added
                        for line in block:
                            if rest_is_code:
                                artifacts_file.write(line)
                            else:
                                try:

                                    if line.startswith("title"): # Write the title to the .md and start reading in
                                        parts = line.split(":")
                                        md_file.write("# " + parts[1])
                                        started = True

                                    elif line.startswith("subtitle"): # Write a subtitle to the .md
                                        parts = line.split(":")
                                        md_file.write("# " + parts[1])

                                    elif line.startswith("image"): # Put an image tag in .md, and code in the python
                                        parts = line.split(":")
                                        image_filename = parts[1]
                                        md_file.write(f"\n![]({image_filename})\n\n")
                                        rest_is_code = True

                                    elif line.startswith("skip"): # Skips a certain number of lines
                                        parts = line.split(":")
                                        if len(parts) <= 1:
                                            skip_lines = 1
                                        else:
                                            skip_lines = int(parts[1])

                                    elif line.startswith("reproduces"):
                                        parts = line.split(":")
                                        md_file.write(spinw_tutorial_text % parts[1].strip())

                                    elif line.startswith("capture-stdout"):
                                        # Five hashes and a colon, all of this is hacky
                                        artifacts_file.write(f'print("#####:{stdout_capture_index}")\n')

                                    elif line.startswith("capture-end"):
                                        md_file.write("\n" + md_placeholder + str(stdout_capture_index) + "\n")
                                        # Six hashes and a colon
                                        artifacts_file.write(f'print("######:{stdout_capture_index}")\n')

                                    else:
                                        raise ValueError(f"Do not know what to do with instruction: {line}")

                                except Exception as e:
                                    raise ValueError(f"Could not parse line {line.strip()}") from e
                    case _:
                        raise ValueError(f"Unknown level {level}")

            artifacts_file.write('\n'
                                 'with open("stdout_data.txt", "w") as file:\n'
                                 '    file.write(buffer.getvalue())\n'
                                 'sys.stdout = old_stdout\n\n')

    if not started:
        print("WARNING: '## title:' was not found, output never started")

    # Execute the artefact files, generating images etc, and insert any data from stdout that needs to be done
    old_dir = os.getcwd()
    os.chdir(output_target_dir)
    runpy.run_path("artifacts.py", run_name="__main__")
    os.chdir(old_dir)

    # Do replacement
    ## Get blocks stdout record
    stdout_blocks = defaultdict(list)
    index = -1
    with open(output_target_dir / "stdout_data.txt", 'r') as file:

        for line in file:
            if line.startswith("#####:"):
                parts = line.split(":")
                index = int(parts[1])
            elif line.startswith("######:"):
                index = -1
            else:
                if index != -1:
                    stdout_blocks[index].append(line)


    # do the replacement
    oldlines = []
    with open(output_target_dir / "tutorial.md", 'r') as md_file:
        oldlines = md_file.readlines()


    with open(output_target_dir / "tutorial.md", 'w') as md_file:
        for line in oldlines:
            if line.startswith(md_placeholder):
                block_index = int(line[len(md_placeholder):].strip())
                md_file.write("```text\n")
                for replace_line in stdout_blocks[block_index]:
                    md_file.write(replace_line)
                md_file.write("```\n\n")
            else:
                md_file.write(line)

    break