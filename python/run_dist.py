import sys
import shlex
import pysam
import subprocess
import multiprocessing
import os


def is_nonempty_file(path):
    return os.path.isfile(path) and os.path.getsize(path)


def submit_mdist(tup):
    ss, ks, sketchfns, ofn = tup
    subprocess.check_call("mash dist -p 1 %s > %s" %
                          (" ".join(sketchfns), ofn), shell=True)


def submit_mash_sketch(tup):
    fa, fn, ss, ks = tup
    assert ss < 32
    ss = 1 << ss
    if not is_nonempty_file(fn + ".msh"):
        assert not os.path.isfile(fn + ".msh")
        print("%s does not exist" % fn, file=sys.stderr)
        cstr = "mash sketch -s {ss} -k {ks} -o {fn} {fa}".format(
            **locals())
        print("About to call '%s'." % cstr, file=sys.stderr)
        subprocess.check_call(shlex.split(cstr))
        print("Successfully called command for fn %s" % fn, file=sys.stderr)


ESTIMS = ["original", "improved", "ertl_ml", "ertl_joint"]
ESTIM_DICT = dict(zip(ESTIMS, ["-E", "-I", "-J", ""]))


def submit_distcmp_call(tup):
    k, n, opath, paths, redo, estim_method = tup
    estim_flag = ESTIM_DICT[estim_method]
    if any(not os.path.isfile(path) for path in paths):
        raise Exception("The files are NOT in the computer! %s" %
                        ", ".join(paths))
        if not redo:
            if is_nonempty_file(opath):
                print("%s has run and redo is set to false. "
                      "Continuing" % opath, file=sys.stderr)
                return
            else:
                print("%s has run, but the file is empty."
                      " Redoing!" % opath, file=sys.stderr)
    cstr = "distcmp %s -So%s -s%s -mn%i -k%i %s" % (
        estim_flag, opath, opath + ".sum", n, k, ' '.join(paths))
    print("Calling '%s'" % cstr, file=sys.stderr)
    subprocess.check_call(shlex.split(cstr))


def make_sketch_fn(fn):
    return fn + ".mash"


def x31_hash(x):
    ret = ord(x[0])
    for i in x[1:]:
        ret = (ret << 5) - ret + ord(i)
    ret = ret & 18446744073709551615  # (1 << 64) - 1
    return ret


def make_hash(x):
    from functools import reduce
    from operator import xor
    print("paths = %s" % x, file=sys.stderr)
    return reduce(xor, map(x31_hash, x))


def makefn(x, y, z, mashify):
    hashval = make_hash(x)
    if not mashify:
        return "experiment_%i_%x_genomes.k%i.n%i.out" % (len(x), hashval, y, z)
    return "mashed_experiment_%i_%x_genomes.k%i.n%i." % (
        len(x), hashval,  y, z)


def main():
    sketch_range = range(10, 24, 1)
    import argparse
    argv = sys.argv
    # Handle args
    superparser = argparse.ArgumentParser(
        description=("This calculates all pairwise distances between "
                     "genomes for all  combinations of parameters."
                     "This does take a while."))
    sp = superparser.add_subparsers()
    shell_parser = sp.add_parser("sketch")
    shell_parser.add_argument("--no-redo", "-n", action="store_true")
    shell_parser.add_argument("--threads", "-p",
                              default=multiprocessing.cpu_count(), type=int)
    shell_parser.add_argument("--use-mash", "-M", action="store_true",
                              help=("Use Mash to calculate distances "
                                    "rather than 'bonsai dist'"))
    shell_parser.add_argument('genomes', metavar='paths', type=str, nargs='+',
                              help=('paths to genomes or a path to a file'
                                    ' with one genome per line.'))
    shell_parser.add_argument("--range-start", default=24, type=int)
    shell_parser.add_argument("--range-end", default=32, type=int)
    py_parser = sp.add_parser(
        "exact",
        description=("Calculates distances natively in Python "
                     "slowly but with obviously no errors"))
    py_parser.add_argument('genomes', metavar='paths', type=str, nargs='+',
                           help=('paths to genomes or a path to a file'
                                 ' with one genome per line.'))
    py_parser.add_argument("--threads", "-p",
                           default=multiprocessing.cpu_count(), type=int)
    py_parser.add_argument("--range-start", default=24, type=int)
    py_parser.add_argument("--range-end", default=32, type=int)
    py_parser.add_argument("--outfile", "-o", default="-")
    if not sys.argv[1:]: sys.argv.append("-h")
    args = superparser.parse_args()
    if argv[1] == "sketch":
        return sketch_main(args)
    elif argv[1] == "exact":
        return exact_main(args)
    else:
       raise Exception("Subcommand required. sketch or exact supported, for sketching or exact calculation.")


def jaccard_index(set1, set2):
    return len(set1 & set2) / float(len(set1 | set2))


def build_kmer_set(tup):
    ks, path = tup
    return {read.sequence[i:i + ks] for read in pysam.FastxFile(path)
            for i in range(len(read.sequence) - ks + 1) if
            all(nuc in 'acgtACGT' for nuc in read.sequence[i:i + ks])}


def exact_main(args):
    kmer_range = range(args.range_start, args.range_end + 1)
    threads = args.threads
    paths = genomes = args.genomes
    if len(genomes) == 1 and os.path.isfile(next(open(genomes[0])).strip()):
        paths = [i.strip() for i in open(genomes[0])]
    Spooooool = multiprocessing.Pool(threads)  # With apologies to Beckett.
    ofp = sys.stdout if args.outfile in ["-", "stdout", "/dev/stdout"] \
        else open(args.outfile, "w")
    ofw = ofp.write  # Cache to reduce the number of lookups

    def pair2tup(path1, path2, *, ks):
        return tuple(sorted([path1, path2]) + [ks])
    ofw("#Path1\tPath2\tKmer Size\tExact Jaccard\n")

    for ks in kmer_range:
        genome_sets = Spooooool.map(build_kmer_set,
                                    ((ks, path) for path in paths))
        kdict = {}
        for i in range(len(genome_sets)):
            for j in range(i+1, len(genome_sets)):
                kdict[tuple(sorted((paths[i], paths[j])) + [ks])] = \
                    jaccard_index(genome_sets[i], genome_sets[j])
        assert len(kdict) == (len(genome_sets) * (len(genome_sets) - 1)) >> 1
        set(not ofw("%s\t%s\t%i\t%f\n" % (*k, v)) for k, v in kdict.items())



def sketch_main(args):
    sketch_range = range(10, 24, 1)
    kmer_range = range(args.range_start, args.range_end + 1)
    threads = args.threads
    redo = not args.no_redo
    paths = genomes = args.genomes
    if len(genomes) == 1 and os.path.isfile(next(open(genomes[0])).strip()):
        paths = [i.strip() for i in open(genomes[0])]
    mashify = args.use_mash
    if any(not os.path.isfile(path) for path in paths):
        raise Exception("The files are NOT in the computer: %s" %
                        ' '.join(path for path in paths if
                                 os.path.isfile(path)))
    Spooooool = multiprocessing.Pool(threads)  # With apologies to Beckett.
    if mashify:
        for ss in sketch_range:
            for ks in kmer_range:
                sketchfns = list(map(lambda x: "%s.k%i.n%i" %
                                     (x[0], x[1], x[2]), ((path, ks, ss) for
                                                          path in paths))
                                 )
                gen = ((path, sketch, ss, ks) for path, sketch in
                       zip(paths, sketchfns))
                while True:
                    # DO IT
                    #
                    # JUST
                    #
                    # DO IT
                    try:
                        Spooooool.map(submit_mash_sketch, gen)
                        break
                    except BlockingIOError:
                        pass
                sketchfns = [i + ".msh" for i in sketchfns]
                fn = makefn(paths, ks, ss, mashify)
                print("Now about to submit dist comparisons", file=sys.stderr)
                todo = []
                thisfns = []
                for i in range(len(paths) - 1):
                    identifier = os.path.basename(paths[i]).split(".")[0]
                    thisfn = "%s.%s.out" % (fn, identifier)
                    if not redo and is_nonempty_file(thisfn):
                        print("Nonempty file, skipping %s" % thisfn,
                              file=sys.stderr)
                        continue
                    else:
                        print("Missing file %s. Generating." % thisfn,
                              file=sys.stderr)
                        todo.append(i)
                        thisfns.append(thisfn)
                gen = ((ss, ks, sketchfns[i:], thisfn) for
                       i, thisfn in zip(todo, thisfns))
                while True:
                    # IF YOU'RE TIRED OF STARTING OVER
                    #
                    # STOP
                    #
                    # GIVING
                    #
                    # UP
                    try:
                        Spooooool.map(submit_mdist, gen)
                        break
                    except BlockingIOError:
                        pass
    else:
        submission_sets = ((ks, ss, makefn(paths, ks, ss, mashify),
                            paths, redo, estim)
                           for ss in sketch_range for ks in kmer_range
                           for estim in ESTIMS)
        Spooooool.map(submit_distcmp_call, submission_sets)
    return 0


if __name__ == "__main__":
    sys.exit(main())
