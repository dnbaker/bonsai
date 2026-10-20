#!/usr/bin/env bash
# CLI smoke test for `bonsai classify` modes: builds khash and compact (-Z)
# databases from two synthetic genomes with a ranked nodes.dmp, then checks
#  - default output is identical with and without -U (classification-only mode);
#  - -U writes a KrakenUniq-style report with the expected header and taxa;
#    -B also writes one b-bit MinHash sketch per clade (and needs -U);
#  - -E adds I: (estimated identity) fields;
#  - -I 0.95,0.85 places 8%-diverged reads of species 562 at its family (543)
#    and leaves error-free reads at 562; khash and compact agree;
#  - -I with a taxonomy file lacking ranks, or a malformed threshold, exits 1
#    with a message (not an abort).
# Usage: test/cli_classify_smoke.sh [path/to/bonsai]
set -euo pipefail
BNS=${1:-./bonsai}
T=$(mktemp -d /tmp/bns_cli_XXXXXX)
trap 'rm -rf "$T"' EXIT
fail() { echo "FAIL: $*" >&2; exit 1; }

# Deterministic random genomes and reads (awk's srand/rand: no external tools)
awk -v seed=11 -v n=30000 -v name=gA 'BEGIN{srand(seed); split("A C G T",b," "); printf(">%s\n", name);
    for(i = 0; i < n; ++i) {printf("%s", b[int(rand() * 4) + 1]); if(i % 80 == 79) printf("\n")} printf("\n")}' > "$T/gA.fa"
awk -v seed=12 -v n=30000 -v name=gB 'BEGIN{srand(seed); split("A C G T",b," "); printf(">%s\n", name);
    for(i = 0; i < n; ++i) {printf("%s", b[int(rand() * 4) + 1]); if(i % 80 == 79) printf("\n")} printf("\n")}' > "$T/gB.fa"
printf '1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n100\t|\t2\t|\torder\t|\n543\t|\t100\t|\tfamily\t|\n561\t|\t543\t|\tgenus\t|\n562\t|\t561\t|\tspecies\t|\n1903414\t|\t100\t|\tfamily\t|\n583\t|\t1903414\t|\tgenus\t|\n584\t|\t583\t|\tspecies\t|\n' > "$T/nodes.dmp"
awk -F'\t' '{print $1 "|" $3}' "$T/nodes.dmp" | sed 's/[[:space:]]//g' > "$T/norank.tax"
printf 'gA\t562\ngB\t584\n' > "$T/seq2tax"
# reads: 6 error-free 5 kb reads of gA (r*), 6 of gA at 8% substitutions (m*), Q40
awk -v seed=13 'BEGIN{srand(seed); split("A C G T",b," ")}
    /^>/ {next} {g = g $0}
    END {
        q = sprintf("%5000s", ""); gsub(/ /, "I", q);
        for(r = 0; r < 6; ++r) {o = int(rand() * (length(g) - 5000)); printf("@r%d\n%s\n+\n%s\n", r, substr(g, o + 1, 5000), q)}
        for(r = 0; r < 6; ++r) {
            o = int(rand() * (length(g) - 5000)); s = substr(g, o + 1, 5000); t = "";
            for(i = 1; i <= 5000; ++i) {c = substr(s, i, 1); if(rand() < 0.08) {do {d = b[int(rand() * 4) + 1]} while(d == c); c = d} t = t c}
            printf("@m%d\n%s\n+\n%s\n", r, t, q)
        }
    }' "$T/gA.fa" > "$T/reads.fq"

"$BNS" build -L -k 31 -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/db.khash" "$T/gA.fa" "$T/gB.fa" > "$T/build.log" 2>&1 || fail "build (khash)"
"$BNS" build -Z -k 31 -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/db.compact" "$T/gA.fa" "$T/gB.fa" > "$T/buildz.log" 2>&1 || fail "build (compact)"
"$BNS" build -k 31 -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/db.default" "$T/gA.fa" "$T/gB.fa" > "$T/buildd.log" 2>&1 || fail "build (default)"
cmp -s "$T/db.default" "$T/db.compact" || fail "build: the default is not the compact (-Z) database"
head -c 6 "$T/db.khash" | grep -q BNSCKM && fail "build -L wrote a compact database"
"$BNS" build -k 31 -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/db.gz" "$T/gA.fa" "$T/gB.fa" > /dev/null 2>&1 && fail "build accepted a gzipped compact database"
"$BNS" build -L -k 31 -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/dbl.gz" "$T/gA.fa" "$T/gB.fa" > /dev/null 2>&1 || fail "build -L .gz"
# feature-count (-f) and taxonomic-depth (-t) minimized builds: prebuild writes
# the phase-1 map, build reads it (<phase1> <out> <genomes>) and leaves it intact
for m in f t; do
    "$BNS" prebuild -$m -k 31 "$T/seq2tax" "$T/nodes.dmp" "$T/p1.$m" "$T/gA.fa" "$T/gB.fa" > "$T/p1.$m.log" 2>&1 || fail "prebuild -$m"
    cp "$T/p1.$m" "$T/p1.$m.orig"
    "$BNS" build -$m -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/p1.$m" "$T/db.min$m" "$T/gA.fa" "$T/gB.fa" > "$T/min$m.log" 2>&1 \
        || { tail -5 "$T/min$m.log"; fail "build -$m"; }
    cmp -s "$T/p1.$m" "$T/p1.$m.orig" || fail "build -$m overwrote its phase-1 input"
    [ -s "$T/db.min$m" ] || fail "build -$m wrote no database"
    head -c 6 "$T/db.min$m" | grep -q BNSCKM && fail "build -$m wrote a compact database"
    "$BNS" classify -k "$T/db.min$m" "$T/nodes.dmp" "$T/reads.fq" > "$T/min$m.out" 2>/dev/null || fail "classify (-$m database)"
    [ "$(awk '$2 ~ /^r/ && $3 == 562' "$T/min$m.out" | wc -l)" -eq 6 ] || fail "-$m database: error-free reads not at 562: $(head -3 "$T/min$m.out")"
done
# Minimizers chosen by the phase-1 score: gC is gA with a SNP every 45 bp under
# another genus (623), so a 60 bp window mixes k-mers shared with gA (LCA: family
# 543) and genome-specific ones. Lex selection keeps some shared ones; -f (fewest
# genomes) and -t (deepest taxon) avoid them. At w > k, density matches lex
# (score ties broken by k-mer hash).
awk 'BEGIN{split("A C G T",b," ")} /^>/ {print ">gC"; next} {g = g $0}
    END {s = ""; for(i = 1; i <= length(g); ++i) {c = substr(g, i, 1); if(i % 45 == 0) for(j = 1; j <= 4; ++j) if(b[j] != c) {c = b[j]; break}; s = s c}
         for(i = 1; i <= length(s); i += 80) print substr(s, i, 80)}' "$T/gA.fa" > "$T/gC.fa"
printf 'gA\t562\ngB\t584\ngC\t623\n' > "$T/seq2tax3"
{ cat "$T/nodes.dmp"; printf '620\t|\t543\t|\tgenus\t|\n623\t|\t620\t|\tspecies\t|\n'; } > "$T/nodes3.dmp"
"$BNS" build -L -k 31 -w 60 -T "$T/nodes3.dmp" -M "$T/seq2tax3" "$T/db3.lex" "$T/gA.fa" "$T/gB.fa" "$T/gC.fa" > /dev/null 2>&1 || fail "build -L -w 60"
shared() { "$BNS" hist "$1" 2>/dev/null | awk -F'\t' '$1 == 543 {n = $2} END {print n + 0}'; }
total() { "$BNS" hist "$1" 2>/dev/null | awk -F'\t' 'NR > 1 {n += $2} END {print n + 0}'; }
lex_shared=$(shared "$T/db3.lex"); lex_total=$(total "$T/db3.lex")
[ "$lex_shared" -gt 100 ] || fail "lex -w 60 kept only $lex_shared shared k-mers (test needs some)"
for m in f t; do
    "$BNS" prebuild -$m -k 31 "$T/seq2tax3" "$T/nodes3.dmp" "$T/p3.$m" "$T/gA.fa" "$T/gB.fa" "$T/gC.fa" > /dev/null 2>&1 || fail "prebuild -$m (3 genomes)"
    "$BNS" build -$m -w 60 -T "$T/nodes3.dmp" -M "$T/seq2tax3" "$T/p3.$m" "$T/db3.$m" "$T/gA.fa" "$T/gB.fa" "$T/gC.fa" > /dev/null 2>&1 || fail "build -$m -w 60"
    n=$(shared "$T/db3.$m"); tot=$(total "$T/db3.$m")
    [ "$n" -lt $((lex_shared / 4)) ] || fail "-$m kept $n shared k-mers (lex: $lex_shared): phase-1 score not used"
    [ "$tot" -lt $((lex_total * 13 / 10)) ] || fail "-$m -w 60 kept $tot k-mers (lex: $lex_total): not minimizing"
done

for db in khash compact; do
    "$BNS" classify -k "$T/db.$db" "$T/nodes.dmp" "$T/reads.fq" > "$T/plain.$db" 2>/dev/null || fail "classify plain ($db)"
    "$BNS" classify -k -U "$T/uniq.$db.tsv" "$T/db.$db" "$T/nodes.dmp" "$T/reads.fq" > "$T/uniq.$db" 2>/dev/null || fail "classify -U ($db)"
    cmp -s "$T/plain.$db" "$T/uniq.$db" || fail "-U changed classification output ($db)"
    head -1 "$T/uniq.$db.tsv" | grep -q '^#taxid	parent	rank	pct_reads	clade_reads' || fail "-U report header ($db)"
    awk -F'\t' '$1 == 562 && $3 == "species" && $5 == 12 {ok = 1} END {exit !ok}' "$T/uniq.$db.tsv" || fail "-U report: 12 reads at species 562 ($db)"
    [ "$(awk '$3 != 562' "$T/plain.$db" | wc -l)" -eq 0 ] || fail "plain calls are not all 562 ($db)"
    mkdir -p "$T/sk.$db"
    "$BNS" classify -k -U "$T/bbmh.$db.tsv" -B "$T/sk.$db" "$T/db.$db" "$T/nodes.dmp" "$T/reads.fq" > "$T/bbmh.$db" 2>/dev/null || fail "classify -U -B ($db)"
    cmp -s "$T/plain.$db" "$T/bbmh.$db" || fail "-B changed classification output ($db)"
    for t in 1 2 100 543 561 562; do [ -s "$T/sk.$db/$t.bbmh" ] || fail "-B: no sketch for clade $t ($db)"; done
    "$BNS" classify -k -E "$T/db.$db" "$T/nodes.dmp" "$T/reads.fq" 2>/dev/null | grep -c 'I:' | grep -qx 12 || fail "-E: I: on every read ($db)"
    "$BNS" classify -k -I 0.95,0.85 "$T/db.$db" "$T/nodes.dmp" "$T/reads.fq" > "$T/placed.$db" 2>/dev/null || fail "classify -I ($db)"
    [ "$(awk '$2 ~ /^r/ && $3 == 562' "$T/placed.$db" | wc -l)" -eq 6 ] || fail "-I: error-free reads not at 562 ($db)"
    [ "$(awk '$2 ~ /^m/ && $3 == 543' "$T/placed.$db" | wc -l)" -eq 6 ] || fail "-I: 8% reads not placed at family 543 ($db)"
done
cmp -s "$T/placed.khash" "$T/placed.compact" || fail "-I: khash and compact differ"
set +e
"$BNS" classify -k -I 0.95,0.85 "$T/db.khash" "$T/norank.tax" "$T/reads.fq" > /dev/null 2> "$T/norank.err"
rc=$?
set -e
[ "$rc" -eq 1 ] || fail "-I with a taxonomy without ranks: exit $rc, want 1"
grep -q 'needs ranks' "$T/norank.err" || fail "-I without ranks: no error message"
"$BNS" classify -k -I bogus "$T/db.khash" "$T/nodes.dmp" "$T/reads.fq" > /dev/null 2>&1 && fail "-I accepted a malformed threshold"
"$BNS" classify -k -B "$T/sk.khash" "$T/db.khash" "$T/nodes.dmp" "$T/reads.fq" > /dev/null 2>&1 && fail "-B accepted without -U"
echo "cli_classify_smoke: all checks passed"
