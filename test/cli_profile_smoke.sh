#!/usr/bin/env bash
# classify -P and bonsai profile: the report recomputed from a profile equals
# classify -U's, profiles of read subsets merge to the whole sample's profile
# byte for byte, and compare/table separate samples of different genomes.
set -euo pipefail
BNS=${1:-./bonsai}
T=$(mktemp -d /tmp/bns_prof_XXXXXX)
[ -n "${KEEP:-}" ] || trap 'rm -rf "$T"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

genome() {
    awk -v seed=$1 -v n=60000 -v name=$2 'BEGIN{srand(seed); split("A C G T",b," "); printf(">%s\n", name);
        for(i = 0; i < n; ++i) {printf("%s", b[int(rand() * 4) + 1]); if(i % 80 == 79) printf("\n")} printf("\n")}'
}
# n reads of length len from a FASTA, seeded
reads() {
    awk -v seed=$2 -v n=$3 -v len=$4 -v pre=$5 'BEGIN{srand(seed)} /^>/ {next} {g = g $0}
        END {q = sprintf("%" len "s", ""); gsub(/ /, "I", q);
             for(r = 0; r < n; ++r) {o = int(rand() * (length(g) - len)); printf("@%s%d\n%s\n+\n%s\n", pre, r, substr(g, o + 1, len), q)}}' "$1"
}
genome 21 gA > "$T/gA.fa"
genome 22 gB > "$T/gB.fa"
printf '1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n100\t|\t2\t|\torder\t|\n543\t|\t100\t|\tfamily\t|\n561\t|\t543\t|\tgenus\t|\n562\t|\t561\t|\tspecies\t|\n1903414\t|\t100\t|\tfamily\t|\n583\t|\t1903414\t|\tgenus\t|\n584\t|\t583\t|\tspecies\t|\n' > "$T/nodes.dmp"
printf 'gA\t562\ngB\t584\n' > "$T/seq2tax"
reads "$T/gA.fa" 31 400 150 a > "$T/a.fq"
reads "$T/gB.fa" 32 400 150 b > "$T/b.fq"
reads "$T/gA.fa" 33 300 150 c > "$T/c.fq"
reads "$T/gB.fa" 34 100 150 d >> "$T/c.fq"      # c: 75% gA, 25% gB
head -n 800 "$T/a.fq" > "$T/a1.fq"; tail -n +801 "$T/a.fq" > "$T/a2.fq"

"$BNS" build -k 31 -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/db" "$T/gA.fa" "$T/gB.fa" > "$T/build.log" 2>&1 || fail build

for mode in hll mh; do
    flag=; [ $mode = mh ] && flag=-m
    "$BNS" classify -k $flag -U "$T/a.$mode.tsv" -P "$T/a.$mode.bprof" "$T/db" "$T/nodes.dmp" "$T/a.fq" > "$T/a.$mode.out" 2>/dev/null || fail "classify -P ($mode)"
    "$BNS" classify -k "$T/db" "$T/nodes.dmp" "$T/a.fq" > "$T/a.plain" 2>/dev/null
    cmp -s "$T/a.plain" "$T/a.$mode.out" || fail "-P changed classification output ($mode)"
    "$BNS" profile report "$T/a.$mode.bprof" "$T/a.$mode.rep" || fail "profile report ($mode)"
    cmp -s "$T/a.$mode.tsv" "$T/a.$mode.rep" || { diff "$T/a.$mode.tsv" "$T/a.$mode.rep" | head; fail "profile report != classify -U ($mode)"; }
    # merge of the two halves == the whole
    for h in a1 a2; do
        "$BNS" classify -k $flag -P "$T/$h.$mode.bprof" "$T/db" "$T/nodes.dmp" "$T/$h.fq" > /dev/null 2>&1 || fail "classify $h ($mode)"
    done
    "$BNS" profile merge -n a.$mode "$T/merged.$mode.bprof" "$T/a2.$mode.bprof" "$T/a1.$mode.bprof" || fail "merge ($mode)"
    cmp -s "$T/merged.$mode.bprof" "$T/a.$mode.bprof" || fail "merge(halves) != whole ($mode)"
    for s in b c; do
        "$BNS" classify -k $flag -P "$T/$s.$mode.bprof" "$T/db" "$T/nodes.dmp" "$T/$s.fq" > /dev/null 2>&1 || fail "classify $s ($mode)"
    done
done
# -B from a MinHash profile: same sketches as classify -U -B
mkdir -p "$T/sk1" "$T/sk2"
"$BNS" classify -k -U "$T/x.tsv" -B "$T/sk1" -P "$T/x.bprof" "$T/db" "$T/nodes.dmp" "$T/a.fq" > /dev/null 2>&1 || fail "classify -B -P"
"$BNS" profile report -B "$T/sk2" "$T/x.bprof" > /dev/null || fail "profile report -B"
for f in "$T"/sk1/*.bbmh; do cmp -s "$f" "$T/sk2/$(basename "$f")" || fail "profile -B sketch $(basename "$f") differs"; done
"$BNS" profile report -B "$T/sk2" "$T/a.hll.bprof" > /dev/null 2>&1 && fail "profile report -B accepted a HyperLogLog profile"

"$BNS" profile compare -K "$T/k.tsv" "$T/a.mh.bprof" "$T/a1.mh.bprof" "$T/b.mh.bprof" "$T/c.mh.bprof" > "$T/cmp.tsv" || fail compare
row() { awk -F'\t' -v a="$1" -v b="$2" '$1 == a && $2 == b' "$T/cmp.tsv"; }
# a vs b: disjoint species (bray-curtis 1); 562 and 584 are 6 unit edges apart
row a.mh b.mh | awk -F'\t' '{exit !($5 == 1 && $6 == 1 && $8 > 5.9 && $8 < 6.01)}' || fail "a vs b distances: $(row a.mh b.mh)"
row a.mh a1.mh | awk -F'\t' '{exit !($5 < 0.02 && $8 < 0.2)}' || fail "a vs a1 distances: $(row a.mh a1.mh)"
row a.mh c.mh | awk -F'\t' '{exit !($5 > 0.15 && $5 < 0.35)}' || fail "a vs c bray-curtis: $(row a.mh c.mh)"
# k-mers of species 562: a1 is half of a's reads (contained in a); a and c sample the same genome
awk -F'\t' '$1 == 562 && $3 == "a.mh" && $4 == "a1.mh" {ok = ($12 > 0.95)} END {exit !ok}' "$T/k.tsv" || fail "a1 k-mers not contained in a: $(awk -F'\t' '$1 == 562' "$T/k.tsv" | head -3)"
awk -F'\t' '$1 == 584 && $3 == "a.mh" && $4 == "b.mh"' "$T/k.tsv" | grep -q . && fail "-K reported 584 shared by a and b"
"$BNS" profile compare -r species "$T/a.hll.bprof" "$T/b.hll.bprof" > /dev/null || fail "compare -r species"
"$BNS" profile compare "$T/a.hll.bprof" "$T/b.mh.bprof" > /dev/null 2>&1 && fail "compare accepted HyperLogLog vs MinHash profiles"
"$BNS" profile table -r species "$T/a.mh.bprof" "$T/b.mh.bprof" "$T/c.mh.bprof" > "$T/tab.tsv" || fail table
awk -F'\t' '$1 == 562 {exit !($3 == 400 && $4 == 0 && $5 >= 290)}' "$T/tab.tsv" || fail "table: $(cat "$T/tab.tsv")"
# index + search: a's nearest sample (besides itself) is a1, then c
"$BNS" profile index "$T/idx.bidx" "$T/a.mh.bprof" "$T/b.mh.bprof" > /dev/null || fail index
"$BNS" profile index -a "$T/idx.bidx" "$T/a1.mh.bprof" "$T/c.mh.bprof" > /dev/null || fail "index -a"
"$BNS" profile search -n 4 "$T/idx.bidx" "$T/a.mh.bprof" > "$T/search.tsv" || fail search
[ "$(awk -F'\t' '$2 == 1 {print $3}' "$T/search.tsv")" = a.mh ] || fail "search: self not first: $(cat "$T/search.tsv")"
[ "$(awk -F'\t' '$2 <= 3 {printf("%s ", $3)}' "$T/search.tsv")" = "a.mh a1.mh c.mh " ] || fail "search order: $(cat "$T/search.tsv")"
"$BNS" profile search -M containment -n 2 "$T/idx.bidx" "$T/a1.mh.bprof" | awk -F'\t' '$2 == 2 {exit !($3 == "a.mh" && $7 > 0.95)}' \
    || fail "search -M containment: a1 not contained in a"
"$BNS" profile index -a "$T/idx.bidx" "$T/a.hll.bprof" > /dev/null 2>&1 && fail "index -a accepted a HyperLogLog profile into a MinHash index"

# diff: group B has less species 562 (abundance); group S is 562 from a 3%
# divergent strain at the same abundance (k-mer content, not reads)
awk -v seed=41 'BEGIN{srand(seed); split("A C G T",b," ")} /^>/ {print ">gA_strain"; next}
    {s = ""; for(i = 1; i <= length($0); ++i) {c = substr($0, i, 1); if(rand() < 0.03) {do {d = b[int(rand() * 4) + 1]} while(d == c); c = d} s = s c} print s}' "$T/gA.fa" > "$T/gA2.fa"
: > "$T/groups.tsv"; : > "$T/strain.tsv"
for i in 1 2 3 4; do
    reads "$T/gA.fa" $((100 + i)) $((300 + 10 * i)) 150 x > "$T/A$i.fq"; reads "$T/gB.fa" $((110 + i)) 100 150 y >> "$T/A$i.fq"
    reads "$T/gA.fa" $((120 + i)) $((100 + 10 * i)) 150 x > "$T/B$i.fq"; reads "$T/gB.fa" $((130 + i)) 300 150 y >> "$T/B$i.fq"
    reads "$T/gA.fa" $((140 + i)) 1500 150 x > "$T/SA$i.fq"
    reads "$T/gA2.fa" $((150 + i)) 1500 150 x > "$T/SB$i.fq"
    for s in A B SA SB; do
        "$BNS" classify -k -m -P "$T/$s$i.bprof" "$T/db" "$T/nodes.dmp" "$T/$s$i.fq" > /dev/null 2>&1 || fail "classify $s$i"
    done
    printf '%s\tA\n' "$T/A$i.bprof" >> "$T/groups.tsv"; printf '%s\tB\n' "$T/B$i.bprof" >> "$T/groups.tsv"
    printf '%s\tA\n' "$T/SA$i.bprof" >> "$T/strain.tsv"; printf '%s\tS\n' "$T/SB$i.bprof" >> "$T/strain.tsv"
done
"$BNS" profile diff -r species "$T/groups.tsv" > "$T/diff.tsv" || fail diff
grep -q '^#permanova	tree_emd	.*relabelings=70' "$T/diff.tsv" || fail "diff: no community PERMANOVA over 70 relabelings: $(head -3 "$T/diff.tsv")"
awk -F'\t' '$1 == 562 {exit !($8 < 0.05 && $7 < -1)}' "$T/diff.tsv" || fail "diff: 562 not lower in B: $(awk -F'\t' '$1 == 562' "$T/diff.tsv")"
awk -F'\t' '$1 == 584 {exit !($8 < 0.05 && $7 > 1)}' "$T/diff.tsv" || fail "diff: 584 not higher in B: $(awk -F'\t' '$1 == 584' "$T/diff.tsv")"
"$BNS" profile diff -K -r species "$T/strain.tsv" > "$T/sdiff.tsv" || fail "diff -K"
awk -F'\t' '$1 == 562 {exit !($8 > 0.2 && $13 < 0.05 && $11 < $10)}' "$T/sdiff.tsv" \
    || fail "diff -K: strain shift not found by k-mers (or reads called it): $(awk -F'\t' '$1 == 562' "$T/sdiff.tsv")"
printf '%s\tA\n%s\tB\n%s\tC\n' "$T/A1.bprof" "$T/B1.bprof" "$T/SA1.bprof" > "$T/three.tsv"
"$BNS" profile diff "$T/three.tsv" > /dev/null 2>&1 && fail "diff accepted three groups"
echo "cli_profile_smoke: all checks passed"
