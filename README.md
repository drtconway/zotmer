# Zotmer
An extensible k-mer analysis workbench written in pure Python.

It runs under the standard Python interpreter, but we recommend using it
with the pypy JIT compiling implementation (http:://pypy.org/).  The pypy
implementation deals well with the blocks of bit-fiddling code which
tend to be in the inner loops, whereas the standard Python is pretty slow.

Installing zotmer
-----------------

Installing under the standard Python environment should be easy:

    $ sudo pip install setuptools
    $ sudo pip install docopt
    $ sudo pip install pykmer
    $ python setup.py build
    $ python setup.py test
    $ sudo python setup.py install
    $ sudo python setup.py install_lib

It's a bit trickier under the Pypy environment because you might
need to install pip under pypy:

    $ wget https://bootstrap.pypa.io/get-pip.py
    $ sudo pypy get-pip.py
    $ sudo pypy -m pip install setuptools
    $ sudo pypy -m pip install docopt
    $ sudo pypy -m pip install pykmer
    $ pypy setup.py build
    $ pypy setup.py test
    $ sudo pypy setup.py install
    $ sudo pypy setup.py install_lib

What is zotmer?
---------------

Zotmer is a collection of k-mer manipulation programs designed to
work together for genomic analysis. Although developed with bacterial
genomics in mind (where individual genomes are generally modest in
size, but where we typically want to compare many genomes, or
contextualize a new genome in terms of those already seen), most
of the methods and the individual programs will work for more complex
genomes.

Currently, the zotmer workbench (accessed through the command _zot_,
provides commands for the following manipulations:
    * k-merizing FASTA/FASTQ files
    * printing out k-mers from a compressed set
    * trimming low frequency k-mers
    * printing the k-mer frequency information
    * computing various genomic distances between sets of k-mers
    * randomly/deterministically sampling k-mers from sets of k-mers
    * merging sets of k-mers
    * finding k-mers likely to denote variants
    * projecting a set of k-mers on to another set
    * scanning k-mer sets to extract alleles based on a set of queries

Why Pure Python?
----------------

Good question. Glad you asked!

Python is a pretty nice language for writing genomics analysis
code. Generators, list support, and so on make it easy to throw together
methods. Combining these things with the computational tricks available
when representing k-mers as integers, makes for a nice expressive
programming environment.

