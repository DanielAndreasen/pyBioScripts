#! /home/Bea/miniconda/bin/python
import argparse


def parser():
    parser = argparse.ArgumentParser(description="Work with .fasta files")
    parser.add_argument("input", help="File name of .fasta file")
    parser.add_argument("-o", "--output", default=None, help="Preffix of output file")
    parser.add_argument("-n", "--nexus", help="Turn off nexus conversion", default=True, action="store_false")
    parser.add_argument("-p", "--phylips", help="Turn off phylips conversion", default=True, action="store_false")
    return parser.parse_args()


def namedit(oldname):
    """
    Change name of key to more readable format
    """
    newname = oldname.replace(">", "")
    if "|" not in newname:
        return newname
    newname = newname.split("|")[3]
    newname = newname.split(".")[0]
    return newname


def pad_missing(a):
    counter = 0
    while a.startswith('-'):
        counter += 1
        a = a[1:]
    if counter:
        a = '?' * counter + a

    a = a[::-1]
    counter = 0
    while a.startswith('-'):
        counter += 1
        a = a[1:]
    if counter:
        a = '?' * counter + a
    return a[::-1]


def phylips(ntax, nchar, matrix, fname):
    """
    Save into .phy format from .fasta format
    """
    output = " {} {}\n".format(ntax, nchar)
    namelength = 0
    for name in matrix.keys():
        namelength = max(namelength, len(name))
    namelength += 1
    for key, value in matrix.items():
        output += '{:nnn}{}\n'.replace('nnn', str(namelength)).format(key, pad_missing(value))

    for suffix in [".fa", ".fas", ".fasta"]:
        if fname.endswith(suffix):
            fname = fname.replace(suffix, ".phy")

    with open(fname, "w") as f:
        f.writelines(output)


def nexus(ntax, nchar, matrix, fname):
    output = "#NEXUS\n"
    output += "begin data;\n"
    output += "\tdimensions ntax={} nchar={};\n".format(ntax, nchar)
    output += "\tformat datatype=dna missing=? gap=-;\n"
    output += "\tmatrix\n"
    for key, value in matrix.items():
        output += "\t{}\t{}\n".format(key, pad_missing(value))
    output += "\n;\nend;\n\n"
    output += "begin mrbayes;\n"
    output += "outgroup HQ224958;\n"
    output += "set autoclose=yes nowarn=yes;\n"
    output += "lset nst=6 rates=gamma ngammacat=4;\n"
    output += "prset applyto=(all) pinvar=fixed(0.4210) statefreqpr=fixed(0.3585, 0.1240, 0.1938, 0.3238) shape=fixed(0.6730) revmatpr=fixed(1.0000, 8.1005, 2.0532, 2.0532, 8.1005, 1.0000);\n"
    output += "mcmcp ngen= 10000000 relburnin=yes burninfrac=0.25 printfreq=1000 samplefreq=1000 nchains=4 savebrlens=yes;\n"
    output += "mcmc;\n"
    output += "sump burnin=2500;\n"
    output += "sumt burnin=2500;\n"

    for suffix in [".fa", ".fas", ".fasta"]:
        if fname.endswith(suffix):
            fname = fname.replace(suffix, ".nex")

    with open(fname, "w") as f:
        f.writelines(output)


args = parser()
fname = args.input
with open(fname, "r") as lines:
    d = {}
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            key = namedit(line)
            value = None
        elif line == "":
            value = None
        else:
            value = line
            nchar = len(value)
        if value is not None:
            d[key] = value
ntax = len(d.keys())

fout = args.output
if fout is None:
    fout = fname
else:
    fout += ".fa"

if args.phylips:
    phylips(ntax, nchar, d, fout)
if args.nexus:
    nexus(ntax, nchar, d, fout)
