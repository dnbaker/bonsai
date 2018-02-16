#ifndef _CLADE_HEADER_H__
#define _CLADE_HEADER_H__
#include <string>
#include <iostream>
#include <cassert>
#include <unordered_map>
#include <cstring>
namespace emp {

enum class ClassLevel:int {
    SUPERKINGDOM     = 2,
    KINGDOM          = 3,
    SUBKINGDOM       = 4,
    SUPERPHYLUM      = 5,
    PHYLUM           = 6,
    SUBPHYLUM        = 7,
    SUPERCLASS       = 8,
    CLASS            = 9,
    SUBCLASS         = 10,
    INFRACLASS       = 11,
    COHORT           = 12,
    SUPERORDER       = 13,
    ORDER            = 14,
    SUBORDER         = 15,
    INFRAORDER       = 16,
    PARVORDER        = 17,
    SUPERFAMILY      = 18,
    FAMILY           = 19,
    SUBFAMILY        = 20,
    TRIBE            = 21,
    SUBTRIBE         = 22,
    GENUS            = 23,
    SUBGENUS         = 24,
    SPECIES_GROUP    = 25,
    SPECIES_SUBGROUP = 26,
    SPECIES          = 27,
    SUBSPECIES       = 28,
    VARIETAS         = 29,
    FORMA            = 30,
    ROOT             = 1,
    NO_RANK          = 0
};

extern const char *classlvl_arr[31];
extern const std::unordered_map<std::string, ClassLevel> classlvl_map;

#define LINE_LVL_OFFSET 0


} // namespace emp

#endif // #ifndef _CLADE_HEADER_H__

