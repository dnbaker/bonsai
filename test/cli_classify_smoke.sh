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

"$BNS" build -k 31 -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/db.khash" "$T/gA.fa" "$T/gB.fa" > "$T/build.log" 2>&1 || fail "build (khash)"
"$BNS" build -Z -k 31 -T "$T/nodes.dmp" -M "$T/seq2tax" "$T/db.compact" "$T/gA.fa" "$T/gB.fa" > "$T/buildz.log" 2>&1 || fail "build (compact)"

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
