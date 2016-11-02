Database Building:
    1. Count the cardinality for windowed sets of genomes for all three scoring schemes.
    2. Build the final jellyfish database hash table.
    3. Provide binary dumps to files for jellyfish databases. (See jellyfish/unit_tests.)
    4. Load binary dumps for jellyfish database, build hash table from the binary dump.
Classification:
    1. Add classification with multiple strategies and supporting multiple databases. (Simple voting, use tree)
    2. Add quantifcation only strategy.
    3. /*later*/Add support for multi-record fastas to work with transcriptome
