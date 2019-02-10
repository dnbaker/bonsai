#ifndef _CLADE_HEADER_H__
#define _CLADE_HEADER_H__
#include <string>
#include <iostream>
#include <cassert>
#include <unordered_map>
#include <cstring>
namespace bns {

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

static const char *classlvl_arr[31] {
    "no rank",
    "root",
    "superkingdom",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "cohort",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "species_group",
    "species_subgroup",
    "species",
    "subspecies",
    "varietas",
    "forma",
};
static const std::unordered_map<std::string, ClassLevel> classlvl_map {
    {"no rank",          ClassLevel::NO_RANK},
    {"root",             ClassLevel::ROOT},
    {"superkingdom",     ClassLevel::SUPERKINGDOM},
    {"kingdom",          ClassLevel::KINGDOM},
    {"subkingdom",       ClassLevel::SUBKINGDOM},
    {"superphylum",      ClassLevel::SUPERPHYLUM},
    {"phylum",           ClassLevel::PHYLUM},
    {"subphylum",        ClassLevel::SUBPHYLUM},
    {"superclass",       ClassLevel::SUPERCLASS},
    {"class",            ClassLevel::CLASS},
    {"subclass",         ClassLevel::SUBCLASS},
    {"infraclass",       ClassLevel::INFRACLASS},
    {"cohort",           ClassLevel::COHORT},
    {"superorder",       ClassLevel::SUPERORDER},
    {"order",            ClassLevel::ORDER},
    {"suborder",         ClassLevel::SUBORDER},
    {"infraorder",       ClassLevel::INFRAORDER},
    {"parvorder",        ClassLevel::PARVORDER},
    {"superfamily",      ClassLevel::SUPERFAMILY},
    {"family",           ClassLevel::FAMILY},
    {"subfamily",        ClassLevel::SUBFAMILY},
    {"tribe",            ClassLevel::TRIBE},
    {"subtribe",         ClassLevel::SUBTRIBE},
    {"genus",            ClassLevel::GENUS},
    {"subgenus",         ClassLevel::SUBGENUS},
    {"species group",    ClassLevel::SPECIES_GROUP},
    {"species subgroup", ClassLevel::SPECIES_SUBGROUP},
    {"species",          ClassLevel::SPECIES},
    {"subspecies",       ClassLevel::SUBSPECIES},
    {"varietas",         ClassLevel::VARIETAS},
    {"forma",            ClassLevel::FORMA},
};

#define LINE_LVL_OFFSET 0


} // namespace bns

#endif // #ifndef _CLADE_HEADER_H__

