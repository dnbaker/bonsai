import random


def gen_seq(nnucs, alphabet="ACGT"):
    return "".join(random.choice(alphabet) for i in range(nnucs))


def gen_fasta(seq, name):
    if not name:
        raise Exception("name must not be null.")
    return ">%s\n%s\n" % (name, seq)


def make_genomes(nunique_per_leaf=13000, nshared_per_subgroup=5000,
                 stuff=0):
    return None


def chunker(iterable, chunk_size):
    return (list(iterable[i:i+chunk_size]) for i in
            range(0, len(iterable), max(1, chunk_size)))


def write_nameid_map(name_id_map, filenames):
    with open(name_id_map, "w") as f:
        for file in filenames:
            name = next(open(file)).split('|')[0][1:]
            taxid = int(name)
            f.write("%s\t%i\n" % (name, taxid))


class SeqId:
    def __init__(self, sequence, id):
        self.seq = sequence
        self.id = id
        assert isinstance(self.seq, str)
        assert isinstance(self.id, int)

        # Keeps track of which taxonomic groups it belongs to.
        self.subsets = {}

        # Keeps track of which untracked taxonomic subgroup it belongs to.
        # This lets makes it easier to reconstruct the taxonomy
        self.subgroups = {}

    def taxid(self):
        # Reserving 0 for null and 1 for the root of the tree.
        return self.id + 2

    def filename(self):
        return "synthseq.%i.fa" % self.taxid()

    def write(self, path):
        def dict2str(d):
            return ';'.join("%s,%s" % (k, v) for k, v in d.items())
        with open(path, "w") as f:
            f.write(">%i|subsets:%s|subgroups:%s|\n" %
                    (self.taxid(),
                     dict2str(self.subsets), dict2str(self.subgroups)))
            for chunk in range(0, len(self.seq), 80):
                f.write(self.seq[chunk:chunk + 80])
                f.write('\n')


__all__ = ["gen_seq", "make_genomes", "gen_fasta", "chunker", "SeqId"]
