import sys
import os

sys.path.insert(0, '../')
from src import *

# Test cases for generating synthetic reference datasets
# Reference 1
#   [26786528]                                                                                                 [26786635]
# chr1: ACCTCTCACTCCTGCCTGGTGTTCCAACCCGTTCTGTGGCCAGAGTATACATTTTGGAACCTCTTCGAGGCCATCCTGCAGTTCCAGATGAACCATAGCGTGCTTCAG
# 1 - base indexing

# Reference 2
#   [116403940]
# chr7: AGTGCAGCTTTTTTATGTCCCTTCAGTAGTTTCTACTTCAAATAATA
# 1 - base indexing

# Substitution of bases
# chr1:26786578-26786603
# REF: ATTTTGGAACCTCTTCGAGGCCATC
# ALT: ATTTTGGAACCGCTTCGAGGCCATC

# chr1    26786590        T    C    SNP

# Deletion of bases
# chr1:26786578-26786603
# REF: ATTTTGGAACCTCTTCGAGGCCATC
# ALT: AGGAACCTCTTCGAGGCCATC

# chr1    26786578    26786582            DEL

# Insertion of bases
# chr1:26786578-26786603
# REF: ATTTTGGAACCTCTTCGAGGCCATC
# ALT: ATTTTGGAACCTCTTCGAGGCCTTAGATC

# chr1    26786600    26786600        TTAG    INS

# Translocation of bases
# chr1: 26786540-26786565
# REF: TGCCTGGTGTTCCAACCCGTTCTGT
# ALT: TGCCTTCCAGATGAACCATAGCG

# chr1: 26786603-26786628
# REF: CTGCAGTTCCAGATGAACCATAGCG
# ALT: CTGCAGTGGTGTTCCAACCCGTTCTGT

# chr1    26786545    chr1    26786610        TRANS
########################################################
# Test loading a variant definition for a substitution #
########################################################
reference = 'ATTTTGGAACCTCTTCGAGGCCATC'
alternate = 'ATTTTGGAACCGCTTCGAGGCCATC'
# Variant Defintion
# chr1    26786590        T    C    SNP

var = SubstitutionVariantDefinition('chr1', '26786590', 'T', 'G', var_type="SNP")
var.set_relative_start(26786578)

assert(var.modify(reference) == alternate)
assert(var.chromosome == 'chr1')
assert(var.position == 26786590)
assert(var.ref == 'T')
assert(var.alt == 'G')

assert(var.alt not in ['C', 'T', 'A'])

####################################################
# Test loading a variant definition for a deletion #
####################################################

reference = 'ATTTTGGAACCTCTTCGAGGCCATC'
alternate = 'AXXXXGGAACCTCTTCGAGGCCATC'

var = dict()

# Variant Definition
# chr1    26786578    26786582            DEL
var['case_3'] = DeletionVariantDefinition('chr1', '26786578', '26786582', var_type="DEL")
var['case_3'].set_relative_start(26786578)

assert(var['case_3'].modify(reference) == alternate)

alternate = 'XXXXXXXXXXXXXXXXXXXXXXXXX'
var['case_1'] = DeletionVariantDefinition('chr1', '26786570', '26786605', var_type="DEL")
var['case_1'].set_relative_start(26786578)

assert(var['case_1'].modify(reference) == alternate)

alternate = 'XXXXXGGAACCTCTTCGAGGCCATC'

var['case_2'] = DeletionVariantDefinition('chr1', '26786577', '26786582', var_type="DEL")
var['case_2'].set_relative_start(26786578)

assert(var['case_2'].modify(reference) == alternate)

######################################################
# Test Loading a variant definition for an insertion #
######################################################

reference = 'ATTTTGGAACCTCTTCGAGGCCATC'
alternate = 'ATT_ATCG_TTGGAACCTCTTCGAGGCCATC'


var = dict()

var['case_3'] = InsertionVariantDefinition('chr1', '26786580', 4, 'ATCG', var_type="INS")
var['case_3'].set_relative_start(26786578)

assert(var['case_3'].modify(reference) == alternate)

var['case_1'] = InsertionVariantDefinition('chr1', '26786577', 5, var_type="INS")
var['case_1'].set_relative_start(26786578)

assert(var['case_1'].modify(reference)[0] == '_' and var['case_1'].modify(reference)[6] == '_')

var['case_2'] = InsertionVariantDefinition('chr1', '26786603', 2, var_type="INS")
var['case_2'].set_relative_start(26786578)

assert(var['case_2'].modify(reference)[len(reference)] == '_' and var['case_2'].modify(reference)[len(reference) + 3] == '_')

var['case_1'].modify(var['case_3'].modify(var['case_2'].modify(reference)))

#########################################################
# Test loading a variant definition for a translocation #
#########################################################

reference1 = 'ATTTTGGAACCTCTTCGAGGCCATC'
reference2 = 'AGTGCAGCTTTTTTATGTCCCTTCA'
reference3 = 'CCTTTATGACGTTGTTGCTGCCTAA'

alternate1 = 'ATTGCTTTTTTATGTCCCTTCA'
alternate2 = 'AGTGCATTGGAACCTCTTCGAGGCCATC'

rel_start_1 = 26786578
rel_start_2 = 116403940
rel_start_3 = 55558958


var = dict()
var['case_1'] = TranslocationVariantDefinition('chr1', '26786580', 'chr7', '116403945', var_type='TRANS')
var['case_1'].set_relative_start(rel_start_1, rel_start_2)

#reference_1_mod, reference_2_mod = var['case_1'].modify(reference1, reference2)
#print reference_1_mod
#print reference_2_mod
#print var['case_1']._apply_mods(
#    *var['case_1']._current_mods(
#        var['case_1'].modify(reference1, reference2)[0]))


var['complex_case_1'] = InsertionVariantDefinition('chr1', '26786603', 2, 'AA', var_type="INS")
var['complex_case_1'].set_relative_start(26786578)


#print var['complex_case_1'].modify(reference1)
unmod, mods = var['case_1']._current_mods(var['complex_case_1'].modify(reference1))
index = var['case_1'].position_1 - var['case_1']._rel_start_1 +1
m = MOD('TRANS',
        index,
        effect='10-TRANS')
mods.append(m)
#print var['case_1']._apply_mods(unmod, mods, tags_only=True)
mod_ref1, mod_ref2 = var['case_1'].modify(var['complex_case_1'].modify(reference1), reference2)

var['case_2'] = TranslocationVariantDefinition('chr1', '26786590', 'chr4', '55558968', var_type='TRANS')
var['case_2'].set_relative_start(rel_start_1, rel_start_3)



#print var['case_2'].modify(var['case_1'].modify(var['complex_case_1'].modify(reference1), reference2)[0], reference3)
unmod, mods = var['case_2']._current_mods(var['case_2'].modify(var['case_1'].modify(var['complex_case_1'].modify(reference1), reference2)[0], reference3)[0])



#print (var['case_2']._apply_mods(unmod, mods))
#print remove_nesting(var['case_2']._apply_mods(unmod, mods))
#print VariantDefinition.TRANSLATION_TABLE


var['case_3'] = TranslocationVariantDefinition('chr4', '55558968', 'chr1', '26786590',  var_type='TRANS')
var['case_3'].set_relative_start(rel_start_3, rel_start_1)
mod_ref3, mod_ref1 = var['case_3'].modify(reference3, reference1)
#mod_ref1_2, mod_ref3_2 = var['case_2'].modify(mod_ref1, mod_ref3)

unmod, mods = var['case_2']._current_mods(mod_ref3)
print mod_ref1
print mod_ref3
print unmod, mods
print var['case_2']._apply_mods(unmod, mods)
print var['case_2']._apply_mods(unmod, mods)
print var['case_2']._apply_mods(unmod, mods)

unmod, mods = var['case_2']._current_mods(mod_ref1)
print unmod, mods
print var['case_2']._apply_mods(unmod, mods)


print VariantDefinition.TRANSLATION_TABLE

######################################################
# Test loading a variant definition for an inversion #
######################################################

reference1 = 'ATTTTGGAACCTCTTCGAGGCCATC'
reference2 = 'AGTGCAGCTTTTTTATGTCCCTTCA'
reference3 = 'CCTTTATGACGTTGTTGCTGCCTAA'

alternate1 = 'ATTGCTTTTTTATGTCCCTTCA'
alternate2 = 'AGTGCATTGGAACCTCTTCGAGGCCATC'

rel_start_1 = 26786578
rel_start_2 = 116403940
rel_start_3 = 55558958


var = dict()
var['case_1'] = TranslocationVariantDefinition('chr1', '26786580', 'chr7', '116403945', var_type='INV')
var['case_1'].set_relative_start(rel_start_1, rel_start_2)

mod_ref1, mod_ref2 = var['case_1'].modify(reference1, reference2)

unmod, mods = var['case_1']._current_mods(mod_ref1)
print mod_ref1
print mod_ref2
print unmod, mods
print var['case_1']._apply_mods(unmod, mods)

unmod, mods = var['case_1']._current_mods(mod_ref2)
print var['case_1']._apply_mods(unmod, mods)

mod_ref1, mod_ref2 = var['case_1'].modify(reference1, reference2)
unmod_ref1, mods = var['case_1']._current_mods(mod_ref1)
print mods

print mod_ref1

var['case_1']._apply_mods(unmod_ref1, mods)

#print var['case_1'].modify(reference)
#print var['case_2'].modify(reference)
#print var['case_3'].modify(var['case_2'].modify(reference))


#print var['case_1'].modify(var['complex_case_1'].modify(reference1), reference2)
print var['case_1'].TRANSLATION_TABLE
print var['case_1'].TRANS_COUNT
print var['case_1'].TRANSLATION_TABLE
var['case_3'] = TranslocationVariantDefinition('chr1', '26786590', 'chr4', '55558968', var_type='TRANS')
var['case_3'].set_relative_start(rel_start_1, rel_start_3)
var['case_3'].modify(reference1, reference3)
print var['case_3'].TRANS_COUNT

VariantDefinition.TRANS_COUNT

#####################################################
# Test loading a variant definition for copy number #
#####################################################

reference1 = 'ATTTTGGAACCTCTTCGAGGCCATC'

rel_start_1 = 26786578

var = dict()

var['case_1'] = CopynumberVariantDefinition('chr1', '26786580', '26786585', '3', var_type='CNV')
var['case_1'].set_relative_start(rel_start_1)
var['case_1'].modify(reference1)

