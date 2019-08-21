### Search and get sequences from Genbank's Nucleotide database
### and save as individual genbank files and a csv with general


from Bio import Entrez, SeqIO
from glob import glob
import argparse
import sys
import os


class Taxonomy:
    def __init__(self, name):
        self.name = name
        self.taxid = self._get_tax_id()
        if self.taxid is not None:
            self.taxdata = self._get_tax_data()
        else:
            self.taxdata = None

    def _get_tax_id(self):
        s = Entrez.esearch(db='taxonomy', term=self.name)
        try:
            return Entrez.read(s)['IdList'][0]
        except:
            return None

    def _get_tax_data(self):
        s = Entrez.read(Entrez.efetch(id=self.taxid, db='taxonomy'))
        return s[0]['LineageEx']

    @property
    def cl(self):
        return self.extract_data('class')

    @property
    def family(self):
        return self.extract_data('family')

    @property
    def order(self):
        return self.extract_data('order')

    @property
    def species(self):
        return self.extract_data('species')

    def extract_data(self, value):
        if self.taxdata is None:
            return ''
        for l in self.taxdata:
            if l['Rank'] == value:
                if l['ScientificName'] is None:
                    return ''
                return l['ScientificName']
        return ''


def _prepare_search(s):
    if not os.path.isfile(s):
        # Simple search string
        return s
    codes = []
    # We are dealing with a mr. bayes tree file
    with open(s, 'r') as lines:
        for i, line in enumerate(lines):
            if i < 4:
                # skip header
                continue
            line = line.strip()
            if line == ';':
                # We have all the codes we are looking for
                break
            if len(line) != 8:
                continue
            pre, post = line[0:2], line[2:]
            if all(w.isupper() for w in pre) and post.isdigit():
                codes.append(line)
    return ','.join(codes)


class read2csv():
    def __init__(self):
        self.data = []

    def addSeq(self, seq):
        assession = seq.name
        organism = seq.annotations['organism']
        organism2 = ' '.join(organism.split(' ')[0:2])
        sequence = str(seq.seq)

        ref = seq.annotations["references"][0]
        authors = ref.authors
        journal = ref.journal
        reftitle = ref.title

        host = ""
        host2 = ""
        isolate = ""
        country = ""
        country2 = ""
        product = ""
        cl = ""
        family = ""
        order = ""
        species = ""

        for feat in seq.features:
            qual = feat.qualifiers
            if "host" in qual.keys():
                host = "|".join(qual["host"])
                host2 = ' '.join(host.split(' ')[0:2])
                tax = Taxonomy(host2)
                cl = tax.cl
                family = tax.family
                order = tax.order
                species = tax.species
            if "isolate" in qual.keys():
                isolate = "|".join(qual["isolate"])
            if "country" in qual.keys():
                country = "|".join(qual["country"])
                country2 = country.split(':')[0]
            if "product" in qual.keys():
                product = "|".join(qual["product"])

        self.data.append([assession, organism, organism2, host, host2, species, family, order, cl,
                          isolate, country, country2, product, authors, reftitle, journal, sequence])

    def write2csv(self, file, sep=","):
        file = open(file, "w")
        file.write(sep.join(["assession", "organism", "organism2", "host", "host2", "species",
                             "family", "order", "cl", "isolate",
                             "country", "country2", "product", "authors", "reftitle",
                             "journal", "sequence"]) + "\n")
        for line in self.data:
            file.write(sep.join(line) + "\n")

        file.close()




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Get sequences from genbank')

    parser.add_argument('--email', dest='email', action='store',
                        default=None, help='Email for the Entrez',
                        type=str)
    parser.add_argument('--search', dest='search', action='store',
                        default=None, help='Search string (e.g. Agama[orgn])',
                        type=str)
    parser.add_argument('--csv', dest='outcsv', action='store',
                        default=None, help='Output csv file for results"',
                        type=str)
    parser.add_argument('--out', dest='outdir', action='store',
                        default=None, help='Output folder for genbank files"',
                        type=str)
    parser.add_argument('--bsize', dest='batchSize', action='store',
                        default=100, help='The batch size for processing between downloads"',
                        type=int)
    parser.add_argument('--clean', action='store_true', default=False,
                        help='Remove all .gb files')
    args = parser.parse_args()


    if args.email == None:
        print("Need an email for Entrez. See help")
        sys.exit(1)

    if args.search == None:
        print("Need a search string. See help")
        sys.exit(1)

    check = [str(int(x)) for x in [args.outcsv != None and args.outdir != None]]
    if check == "00":
        print("Need a csv file or a directory to output results. See help")
        sys.exit(1)
    elif check == "01":
        print("Warning: Genbank file will be saved in %s but no csv will be produced" % args.outdir)
    elif check == "10":
        print("Warning: the file %s will be generate but no Genbank file will be saved" % args.outcsv)

    search = _prepare_search(args.search)

    batch_size = args.batchSize


    Entrez.email = args.email
    Entrez.tool = "SearchAndGetPythonScript"

    # Get results
    with Entrez.esearch(db="nucleotide", term=search, usehistory="y") as handle:
        results = Entrez.read(handle)

    count = int(results["Count"])
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]

    print("%s records found." % (count))

    csv = read2csv()

    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        print("Downloading records from %i to %i" % (start+1, end))

        fhandle = Entrez.efetch(db="nucleotide", retmode="text",
                                rettype = "gb", retstart=start,
                                retmax=batch_size, webenv=webenv,
                                query_key=query_key)

        data = SeqIO.parse(fhandle, "genbank")

        for record in data:
            idl = record.name
            print("-----> processing %s" % (idl))
            csv.addSeq(record)
            if args.outdir != None:
                SeqIO.write(record, os.path.join(args.outdir, idl + ".gb"),
                            "genbank")
        fhandle.close()

    if args.outcsv != None:
        csv.write2csv(args.outcsv)

    if args.clean:
        for f in glob('*gb'):
            os.remove(f)

    print("Done!")
