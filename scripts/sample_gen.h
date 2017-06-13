//.h:
enum class ClassLevel:int {
    SUPERKINGDOM = 0,
    KINGDOM      = 1,
    SUBKINGDOM   = 2,
    SUPERPHYLUM  = 3,
    PHYLUM       = 4,
    SUBPHYLUM    = 5,
    SUPERCLASS   = 6,
    CLASS        = 7,
    SUBCLASS     = 8,
    INFRACLASS   = 9,
    COHORT       = 10,
    SUPERORDER   = 11,
    ORDER        = 12,
    SUBORDER     = 13,
    INFRAORDER   = 14,
    PARVORDER    = 15,
    SUPERFAMILY  = 16,
    FAMILY       = 17,
    SUBFAMILY    = 18,
    TRIBE        = 19,
    SUBTRIBE     = 20,
    GENUS        = 21,
    SUBGENUS     = 22,
    SPECIES      = 23,
    SUBSPECIES   = 24,
    GROUP        = 25,
    SUBGROUP     = 26,
    VARIETAS     = 27,
    FORMA        = 28,
    ROOT         = -1,
    NO_RANK      = -2
};

extern const char *classlvl_arr[31];
extern const std::unordered_map<std::string, ClassLevel> classlvl_map;

