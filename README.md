# Table of Contents
- [Getting started](#getting-started)
- [Installation](#installation)
  - [Obtaining annotarium](#obtaining-annotarium)
  - [Prerequisites](#prerequisites)
  - [Installation suggestions](#installation-suggestions)
- [How to cite](#how-to-cite)

# Getting started
```
# Download this repository [making sure you have prerequisite Python packages]
git clone https://github.com/zkstewart/annotarium.git

# Generate statistics detailing a FASTA file
python /location/of/annotarium/annotarium.py fasta stats \
    -i $FASTA -o fasta_stats.tsv

# Convert GFF3 annotations into FASTA sequences
python /location/of/annotarium/annotarium.py gff3 to fasta \
    -i $GFF3 -f $FASTA -o sequences_prefix

# Relate any GFF3 feature values to other details in the GFF3 file
python /location/of/annotarium/annotarium.py gff3 to tsv \
    -i $GFF3 -o gff3_ids_map.tsv \
    -forEach mRNA -map ID -to Parent start end # etc...

# Reformat a GFF3 into proper GFF3 standard, with sorting
python /location/of/annotarium/annotarium.py gff3 to gff3 \
    -i $GFF3 -o updated_annotation.gff3 
```

# Installation
## Obtaining annotarium
Download annotarium by cloning the repository as below. It is available as a collection of Python scripts, so no further installation or compilation is necessary.

```
git clone https://github.com/zkstewart/annotarium.git
```

## Prerequisites
annotarium is written in Python and requires a modern version 3. It makes use of several other Python packages which include:
- biopython (https://biopython.org/)
- numpy (https://numpy.org/)
- pandas (https://pandas.pydata.org/)
- ncls (https://pypi.org/project/ncls)

annotarium has been developed on Linux and within Windows Subsystem for Linux (WSL), but at this stage is likely to be fully Windows compatible.

## Installation suggestions
If you are interested in running annotarium, you should ideally set up an Anaconda or Miniconda environment containing a recent Python 3 version. All Python packages listed above are available through community channels including conda-forge.

See the [Installing annotarium wiki page](https://github.com/zkstewart/annotarium/wiki/Installing-annotarium) for more information.

# How to use annotarium
On the command line, you can always ask annotarium to provide help information for its different functions by doing:

```
python /location/of/annotarium.py gff3 -h
python /location/of/annotarium.py fasta -h
... etc ...
```

Otherwise, refer to the [annotarium wiki](https://github.com/zkstewart/annotarium/wiki) for more detailed information on the program.

# How to cite
There are no plans to publish any of the code in this repository. However, if you do find any of this useful in your own work, please feel free to link to this repository in your manuscript.
