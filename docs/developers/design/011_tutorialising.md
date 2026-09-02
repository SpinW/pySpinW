# How to write tutorials

## Basic idea

Write a python file, the comments will be turned into a text in markdown.

Beyond this, there are extra functions to deal with images and other aspects of formatting, 
these are specified with double hash, a single space, the command (always lowercase) and, if there are
parameters, a colon followed by the parameter e.g.

```
## title: The title
```

The `title` command is unique as it is required, and its location specifies the start of the output.

## Commands

 - `## title: TITLE` specifiy the title
 - `## reproduces: SPINW_TUTORIAL` adds a notice saying that it reproduces a spinW tutorial
 - `## subtitle: SUBTITLE` add a subtitle section
 - `## skip: N` skip the following `N` lines in the processed output (text or code), often used to 
  supress commands that block, and replace
 - `## image: FILENAME` insert an image at this location, if there are more `## ` lines after, they will be executed
  when building the documents and can be used to create files, e.g. `## plt.savefig(FILENAME)`
 - `## capture-stdout` starts capturing data sent to stdout, it will be inserted as verbatim text when the capture ends
 - `## capture-end` ends the capture