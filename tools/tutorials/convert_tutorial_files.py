""" Generate files for tutorial_inputs that will go on the website """

import os
import shutil
import traceback
from collections import defaultdict
from pathlib import Path
import runpy
import sys
import re

from .jupyter_files import JupyterFile

tools_dir = Path(__file__).parent

tutorial_input_dir = tools_dir / Path("../../examples/tutorials/tutorial_inputs")
tutorial_output_dir = tools_dir / Path("../../examples/tutorials/tutorial_outputs")
image_extensions = ['.png', '.jpg', '.jpeg', '.svg']

RE_IMAGE_BLOCK = re.compile(r'^#\$\$\$\$:(?P<im>[0-9]*)(.*)^#\$\$\$\$\$:(?P=im)', re.DOTALL | re.MULTILINE)
RE_IMAGE_FILE = re.compile(r'[\'\"](.*)(' + '|'.join(image_extensions) + ')')


def _parse_images():
    with open("artifacts.py", "r") as artifacts_file:
        artifacts_py = artifacts_file.read()
    image_blocks = []
    for block in RE_IMAGE_BLOCK.finditer(artifacts_py):
        image_code = block.group(2)
        # Search image code for _first_ output filename with any extension defined above
        filename = RE_IMAGE_FILE.search(image_code)
        if filename is None:
            raise RuntimeError(f'No image filename found in code block {image_code}')
        filename = filename.group(1) + filename.group(2)
        if 'savefig' in image_code:
            savefig_cmd = re.search(r'(\w*)\.savefig\((.*[\'"]' + filename + r'[\'"].*\))', image_code, re.DOTALL | re.MULTILINE)
            savefig_dark_cmd = 'savefig_dark(' + savefig_cmd.group(1) + ',' + savefig_cmd.group(2)
            savefig_dark_cmd = savefig_dark_cmd.replace(filename, filename.replace('.', '_dark.'))
            artifacts_py = artifacts_py.replace(savefig_cmd.group(2), savefig_cmd.group(2) + '\n' + savefig_dark_cmd)
            filename_dark = filename.replace('.', '_dark.')
            image_blocks.append([f'```{{image}} {filename}\n:class: only-light\n```\n',
                                 f'```{{image}} {filename_dark}\n:class: only-dark\n```\n'])
        elif 'snapshot' in image_code:
            snapshot_cmd = re.search(r'(snapshot\(.*[\'"]' + filename + r'[\'"].*\))', image_code, re.DOTALL | re.MULTILINE).group(1)
            snapshotlight_cmd = snapshot_cmd.replace('snapshot', 'snapshot_light')
            snapshotlight_cmd = snapshotlight_cmd.replace(filename, filename.replace('.', '_light.'))
            snapshotdark_cmd = snapshot_cmd.replace('snapshot', 'snapshot_dark')
            artifacts_py = artifacts_py.replace(snapshot_cmd, snapshotdark_cmd + '\n' + snapshotlight_cmd)
            filename_light = filename.replace('.', '_light.')
            image_blocks.append([f'```{{image}} {filename}\n:class: only-dark\n```\n',
                                 f'```{{image}} {filename_light}\n:class: only-light\n```\n'])
        else:
            image_blocks.append([f"\n![]({filename})\n\n"])

    with open("artifacts.py", "w") as artifacts_file:
        artifacts_file.write(artifacts_py)

    return image_blocks


def run(files):

    spinw_tutorial_link = "[Tutorial %d](../../../spinw/tutorials/tutorial%02d)"
    spinw_tutorial_text = f"\n_This tutorial mirrors MATLAB spinW {spinw_tutorial_link}"
    md_placeholder_stdout = "%%%%% stdout "
    md_placeholder_stderr = "%%%%% stderr "
    md_placeholder_image = "%%%%% image "


    # Clean directory
    if os.path.exists(tutorial_output_dir):

        for path in tutorial_output_dir.iterdir():
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()

    else:
        os.mkdir(tutorial_output_dir)

    # Create files

    for md_file in files:

        try:
            if md_file.startswith("_"):
                continue

            if not md_file.endswith(".py"):
                continue

            base_name = md_file.split(".")[0]

            print(base_name)


            output_target_dir = tutorial_output_dir / base_name

            os.mkdir(output_target_dir)

            # Parse file into blocks based on comment indentation

            comment_level = 0
            blocks = []
            current_block = []
            skip_lines = 0
            with open(tutorial_input_dir / md_file, 'r') as input_file:
                for line in input_file:
                    if line.startswith("## "):
                        # Parses skipped lines here so we don't break code blocks
                        if line[2:].strip().startswith("skip"):
                            parts = line.split(":")
                            skip_lines = 2 if len(parts) <= 1 else int(parts[1]) + 1
                            add_line = '"""skip\n'
                        else:
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

                    if skip_lines == 1:
                        current_block.append('skip"""\n')
                        skip_lines = 0
                    if skip_lines > 1:
                        skip_lines -= 1

                blocks.append((comment_level, current_block))

            blocks = [(level, block) for level, block in blocks if not (level == 2 and len(block) == 0)]

            # for level, block in blocks:
            #     print(f"Level {level}:")
            #     for line in block:
            #         print("    ", line.rstrip())

            # Create the corresponding output files
            started = False
            stdout_capture_index = 0
            stderr_capture_index = 0
            image_index = 0

            jupyter = JupyterFile()

            with open(output_target_dir / "index.md", 'w') as md_file:
                with open(output_target_dir / "artifacts.py", 'w') as artifacts_file:
                    artifacts_file.write("import sys\n"
                                         "import io\n"
                                         "import matplotlib.pyplot as plt\n"
                                         "from tools.tutorials.darkmode import *\n"
                                         "old_stdout = sys.stdout\n"
                                         "old_stderr = sys.stderr\n"
                                         "stdout_buffer = io.StringIO()\n"
                                         "stderr_buffer = io.StringIO()\n"
                                         "sys.stdout = stdout_buffer\n"
                                         "sys.stderr = stderr_buffer\n"
                                         )

                    for level, block in blocks:
                        match level:

                            case 0: # It's code
                                if started:
                                    # Clear any blank lines from md
                                    md_lines = [line for line in block if line.strip() != ""]

                                    # Write to .md
                                    if len(md_lines) > 0:
                                        md_file.write("\n```python\n")

                                        for line in md_lines:
                                            md_file.write(line.replace('"""skip\n','').replace('skip"""\n',''))

                                        md_file.write("```\n")

                                    # Write to artefacts file
                                    for line in block:
                                        artifacts_file.write(line)

                                    jupyter.add_code([l.replace('"""skip\n','').replace('skip"""\n','') for l in md_lines])

                                else:
                                    pass #print(f"'{block}' comes before title tag, skipping")


                            case 1: # It's text
                                if started:
                                    for line in block:
                                        md_file.write(line)

                                    jupyter.add_text(block)
                                else:
                                    pass # print(f"'{block}' comes before title tag, skipping")

                            case 2: # It's a control tag
                                rest_is_code = False # Flag to treat the rest of this block as code to be added
                                is_image_block = False
                                for line in block:
                                    if rest_is_code:
                                        artifacts_file.write(line)
                                    else:
                                        try:

                                            if line.startswith("title"): # Write the title to the .md and start reading in
                                                parts = line.split(":")

                                                jupyter.add_text("# " + parts[1])

                                                md_file.write("# " + parts[1])
                                                md_file.write(
                                                    f"\nSource: "
                                                    f"[python (.py)]({base_name}.py) | "
                                                    f"[jupyter (.ipynb)]({base_name}.ipynb) \n\n")

                                                started = True

                                            elif line.startswith("subtitle"): # Write a subtitle to the .md
                                                parts = line.split(":")

                                                jupyter.add_text("## " + parts[1])

                                                md_file.write("## " + parts[1])

                                            elif line.startswith("image"): # Put an image tag in .md, and code in the python
                                                md_file.write("\n" + md_placeholder_image + str(image_index) + "\n")
                                                artifacts_file.write(f'\n#$$$$:{image_index}\n')
                                                rest_is_code = True
                                                is_image_block = True

                                            elif line.startswith("reproduces"):
                                                indices = [int(v) for v in line.split(":")[1].strip().split(",")]
                                                md_file.write(spinw_tutorial_text % (indices[0], indices[0]))
                                                for n in range(1, len(indices)):
                                                    md_file.write(', ' + spinw_tutorial_link % (indices[n], indices[n]))
                                                md_file.write("_\n\n")

                                            elif line.startswith("capture-stdout"):
                                                # Five hashes and a colon, all of this is hacky
                                                stdout_capture_index += 1
                                                artifacts_file.write(f'print("#####:{stdout_capture_index}")\n')

                                            elif line.startswith("end-capture-stdout"):
                                                md_file.write(
                                                    "\n" + md_placeholder_stdout + str(stdout_capture_index) + "\n")
                                                # Six hashes and a colon
                                                artifacts_file.write(f'print("######:{stdout_capture_index}")\n')

                                            elif line.startswith("capture-stderr"):
                                                stderr_capture_index += 1
                                                artifacts_file.write(
                                                    f'print("#####:{stderr_capture_index}", file=sys.stderr)\n')

                                            elif line.startswith("end-capture-stderr"):
                                                md_file.write(
                                                    "\n" + md_placeholder_stderr + str(stderr_capture_index) + "\n")
                                                artifacts_file.write(
                                                    f'print("######:{stderr_capture_index}", file=sys.stderr)\n')

                                            elif line.startswith('"""skip'):
                                                artifacts_file.write(line)

                                            else:
                                                raise ValueError(f"Do not know what to do with instruction: {line}")

                                        except Exception as e:
                                            raise ValueError(f"Could not parse line {line.strip()}") from e
                                if is_image_block:
                                    artifacts_file.write(f'\n#$$$$$:{image_index}\n')
                                    image_index += 1
                            case _:
                                raise ValueError(f"Unknown level {level}")

                    artifacts_file.write('\n'
                                         'with open("stdout_data.txt", "w", encoding="utf-8") as file:\n'
                                         '    file.write(stdout_buffer.getvalue())\n'
                                         'sys.stdout = old_stdout\n\n')

                    artifacts_file.write('\n'
                                         'with open("stderr_data.txt", "w", encoding="utf-8") as file:\n'
                                         '    file.write(stderr_buffer.getvalue())\n'
                                         'sys.stderr = old_stderr\n\n')

            if not started:
                print("WARNING: '## title:' was not found, output never started")

            # Save the jupyter file
            jupyter.write_notebook(output_target_dir / f"{base_name}.ipynb")


            # Make a file without the order two comments
            with open(output_target_dir / f"{base_name}.py", 'w') as file:
                for level, block in blocks:
                    if level == 0:
                        for line in block:
                            file.write(line.replace('"""skip\n','').replace('skip"""\n',''))
                    if level == 1:
                        for line in block:
                            file.write("#" + line)

            # Execute the artefact files, generating images etc, and insert any data from stdout that needs to be done
            old_dir = os.getcwd()
            os.chdir(output_target_dir)

            # File will probably set these to something else, and if there's an error it will be a problem
            old_stdout = sys.stdout
            old_stderr = sys.stderr

            # Parses the artifacts file and change image blocks to make light and dark images and get file name
            image_blocks = _parse_images()

            try:
                runpy.run_path("artifacts.py", run_name="__main__")
            except Exception as e:
                sys.stdout = old_stdout
                sys.stderr = old_stderr
                traceback.print_exc()
                raise e

            # Reset, in case something else has gone wrong
            sys.stdout = old_stdout
            sys.stderr = old_stderr

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

            stderr_blocks = defaultdict(list)
            index = -1
            with open(output_target_dir / "stderr_data.txt", 'r') as file:

                for line in file:
                    if line.startswith("#####:"):
                        parts = line.split(":")
                        index = int(parts[1])
                    elif line.startswith("######:"):
                        index = -1
                    else:
                        if index != -1:
                            stderr_blocks[index].append(line)



            # do the replacement
            oldlines = []
            with open(output_target_dir / "index.md", 'r') as md_file:
                oldlines = md_file.readlines()


            with open(output_target_dir / "index.md", 'w') as md_file:
                for line in oldlines:
                    if line.startswith(md_placeholder_stdout):
                        block_index = int(line[len(md_placeholder_stdout):].strip())
                        md_file.write("```text\n")
                        for replace_line in stdout_blocks[block_index]:
                            md_file.write(replace_line)
                        md_file.write("```\n\n")
                    elif line.startswith(md_placeholder_stderr):
                        block_index = int(line[len(md_placeholder_stderr):].strip())
                        md_file.write("```text\n")
                        for replace_line in stderr_blocks[block_index]:
                            md_file.write(replace_line)
                        md_file.write("```\n\n")
                    elif line.startswith(md_placeholder_image):
                        block_index = int(line[len(md_placeholder_image):].strip())
                        for replace_line in image_blocks[block_index]:
                            md_file.write(replace_line)
                    else:
                        md_file.write(line)


        except Exception as e:
            traceback.print_exc()
            raise ValueError(f"Problem with '{base_name}: {e}'") from e



if __name__ == "__main__":
    files = [file for file in os.listdir(tutorial_input_dir) if not file.startswith("_")]
    run(files[:1])
