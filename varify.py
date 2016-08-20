#!/usr/bin/python
import os
import sys
import subprocess
import argparse
from src import *
import collections
import bisect
import re



FASTA = collections.namedtuple('FASTA',
                               ['chromosome', 'start', 'end', 'sequence'])

FASTA_REF = collections.namedtuple('FASTA_REF',
                                    ['name', 'seq'])



VARIANT = collections.namedtuple('VARIANT',
    ['chromosome', 'start_position', 'end_position', 'modification'])

REGEX = {'coordinates': '(chr[0-9XY]+):g.([0-9]*)_*([0-9]*)',
         'variant': '([ACGT]+>[]AGCT]+)|(dup[AGCT]+)|(del[ACTG0-9]+)|(ins[ACTG0-9]+)'}


def load_variant(data):
    var_type = data[-1]
    var_def = None
    if var_type == 'SNP':
        var_def = SubstitutionVariantDefinition(*data[:4], var_type=var_type)

    elif var_type == 'DEL':
        var_def = DeletionVariantDefinition(*data[:3], var_type=var_type)

    elif var_type == 'INS':
        var_def = InsertionVariantDefinition(data[0],
                                             data[1],
                                             int(data[2]),
                                             data[3],
                                             var_type=var_type)

    elif var_type == 'TRANS' or var_type=='INV':
        var_def = TranslocationVariantDefinition(*data[:4], var_type=var_type)

    elif var_type == 'CNV':
        var_def = CopynumberVariantDefinition(*data[:4], var_type=var_type)

    return var_def


def append_to_variant_list(variant_list, variant):
    if variant.type == 'INV':
        variant_list['TRANS'].append(variant)
    else:
        variant_list[variant.type].append(variant)


def get_variant_definition(fn, vl):
    with open(fn, 'rb') as f:
        l = f.readline().strip()
        while l:
            data = l.split('\t')
            variant = load_variant(data)
            if variant.type == 'INV':
                vl['TRANS'].append(variant)
            else:
                vl[variant.type].append(variant)

            l = f.readline().strip()

    return vl


def get_sequences(fn):
    seqs = dict()

    with open(fn, 'rb') as f:
        l = f.readline().strip()
        header, sequence = (None, None)
        while l:
            if header is None:
                header = l
            else:
                chromosome, coords = header.split(':')
                start, end = map(int, coords.split('-'))
                sequence = l.upper()
                chromosome = chromosome[1:]
                entry = FASTA(chromosome,
                              start, end,
                              sequence)
                try:
                    seqs[chromosome].append(entry)
                except KeyError:
                    seqs[chromosome] = []
                    seqs[chromosome].append(entry)

                header, sequence = (None, None)

            l = f.readline().strip()

    for k in seqs.keys():
        seqs[k].sort(key=lambda x:x.start)
    return seqs


def search_references(a, pos):
    starts = [x.start for x in a]
    index = bisect.bisect_right(starts, pos)
    return index if index == 0 else index - 1


def parse_transvar_entry(l):
    ''' given a line from transvar - parse the output
    '''
    try:
        chromosome, start, end = re.search(REGEX['coordinates'], l).groups()
    except AttributeError:
        return None
    start = int(start)
    if end == '':
        end = start
    else:
        end = int(end)

    variants = re.search(REGEX['variant'], l).groups()

    return VARIANT(chromosome, start, end,
                   variants)


def parse_transvar(transvar_out_fn, vl):
    entries = []
    with open(transvar_out_fn, 'rb') as f:
        l = f.readline().strip()
        new_entry = False
        var_count = 0
        while l:
            l = l.split('\t')[4]
            if l == 'coordinates(gDNA/cDNA/protein)':
                new_entry = True
                if new_entry and var_count == 0:
                    # no suitable translation for variant identified
                    # write message to bed file...
                    var_count = 0
                l = f.readline().strip()
                continue
            if l == '././.':
                l = f.readline().strip()
                continue
            new_entry = False
            variant = parse_transvar_entry(l)
            if variant is None:
                l = f.readline().strip()
                continue
            variant_effects = format_transvar_entry(variant)
            var_count += len(variant_effects)
            for v in variant_effects:
                append_to_variant_list(vl, load_variant(v.split('\t')))
            
            l = f.readline().strip()
            if l == '':
                l = f.readline().strip()

    return entries


def format_transvar_entry(entry):
    bed_entries = []
    for effect in entry.modification:
        if effect is None:
            continue
        if effect[:3] == 'del':
            # HANDLE DELETIONS
            
            del_length = None
            try:
                del_length = int(effect[3:])
            except ValueError:
                del_length = len(effect[3:])

            if entry.start_position == entry.end_position:
                del_start, del_end = entry.start_position - 2, entry.start_position - 1
            else:
                del_start = entry.start_position - 2
                del_end = entry.end_position - 1
            bed_entries.append(
                '\t'.join(
                    [entry.chromosome,
                     str(del_start),
                     str(del_end),
                     str(del_length),
                     'DEL']))
        elif effect[:3] == 'ins':
            # HANDLE INSERTIONS
            ins_length = None
            ins_seq = None
            try:
                ins_length = int(effect[3:])
                ins_seq = ''
            except:
                ins_length = len(effect[3:])
                ins_seq = effect[3:]
            bed_entries.append(
                '\t'.join(
                    [entry.chromosome,
                     str(entry.start_position),
                     str(ins_length),
                     ins_seq,
                     'INS']))

        elif effect[:3] == 'dup':
            # HANDLE DUPLICATIONS
            ins_length = None
            ins_seq = None
            try:
                ins_length = int(effect[3:])
                ins_seq = ''
            except:
                ins_length = len(effect[3:])
                ins_seq = effect[3:]
            bed_entries.append(
                '\t'.join(
                    [entry.chromosome,
                     str(entry.start_position),
                     str(ins_length),
                     ins_seq,
                     'INS']))

        else:
            # HANDLE MULTI-NUCLEOTIDE SNPS
            ref, alt = effect.split('>')
            for i,v in enumerate(ref):
                bed_entries.append(
                    '\t'.join(
                        [entry.chromosome,
                         str(entry.start_position),
                         v,
                         alt[i],
                         'SNP']))

    return bed_entries



def main(args):
    '''
    '''
    vl = {
        'SNP': [],
        'DEL': [],
        'INS': [],
        'TRANS': [],
        'CNV': []
    }
    if args.transvar_variants: 
        parse_transvar(args.transvar_variants, vl)
    if args.variants:
        get_variant_definition(args.variants, vl)

    seqs = get_sequences(args.reference)


    output_ref = []

    for v in vl['SNP']:
        try:
            index = search_references(seqs[v.chromosome], v.position)
        except KeyError:
            # the event is not covered by the bed file
            continue
        prev_entry = seqs[v.chromosome][index]
        v.set_relative_start(prev_entry.start)

        try:
            new_entry = FASTA(prev_entry.chromosome,
                                  prev_entry.start,
                                  prev_entry.end,
                                  v.modify(prev_entry.sequence.upper()))
        except AssertionError:
            import pdb; pdb.set_trace()

        seqs[v.chromosome][index] = new_entry


    for v in vl['DEL']:
        try:
            index_start = search_references(seqs[v.chromosome], v.start_position)
            index_end = search_references(seqs[v.chromosome], v.end_position)
        except KeyError:
            # the event is not covered
            continue
        if index_start == index_end:
            prev_entries = [seqs[v.chromosome][index_start]]
        else:
            prev_entries = seqs[v.chromosome][index_start:index_end]

        for index, pe in enumerate(prev_entries):
            v.set_relative_start(pe.start)

            new_entry = FASTA(pe.chromosome,
                              pe.start,
                              pe.end,
                              v.modify(pe.sequence.upper()))

            seqs[v.chromosome][index_start + index] = new_entry


    for v in vl['INS']:
        try:
            index = search_references(seqs[v.chromosome], v.position)
        except KeyError:
            continue
        prev_entry = seqs[v.chromosome][index]
        v.set_relative_start(prev_entry.start)

        new_entry = FASTA(prev_entry.chromosome,
                          prev_entry.start,
                          prev_entry.end,
                          v.modify(prev_entry.sequence.upper()))
        seqs[v.chromosome][index] = new_entry

    for v in vl['TRANS']:
        try:
            index1 = search_references(seqs[v.chromosome_1], v.position_1)
            index2 = search_references(seqs[v.chromosome_2], v.position_2)
        except KeyError:
            continue

        prev_entry_1 = seqs[v.chromosome_1][index1]
        prev_entry_2 = seqs[v.chromosome_2][index2]

        v.set_relative_start(prev_entry_1.start,
                             prev_entry_2.start)


        mod_ref_1, mod_ref_2 = v.modify(prev_entry_1.sequence.upper(),
                                        prev_entry_2.sequence.upper())

        new_entry_1 = FASTA(prev_entry_1.chromosome,
                            prev_entry_1.start,
                            prev_entry_1.end,
                            mod_ref_1)

        new_entry_2 = FASTA(prev_entry_2.chromosome,
                            prev_entry_2.start,
                            prev_entry_2.end,
                            mod_ref_2)

        seqs[v.chromosome_1][index1] = new_entry_1
        seqs[v.chromosome_2][index2] = new_entry_2

        
    for v in vl['TRANS']:
        try:
            index1 = search_references(seqs[v.chromosome_1], v.position_1)
            index2 = search_references(seqs[v.chromosome_2], v.position_2)
        except KeyError:
            continue

        mod_ref_1 = seqs[v.chromosome_1][index1]
        mod_ref_2 = seqs[v.chromosome_2][index2]

        unmod, mods = v._current_mods(mod_ref_1.sequence)
        mod_ref_1_seqs = set(remove_nesting(v._apply_mods(unmod, mods)))
        counter = 0
        for x in mod_ref_1_seqs:
            counter += 1
            output_ref.append(
                FASTA_REF('{}:{}-{}_TRANS-{}'.format(mod_ref_1.chromosome,
                                                     mod_ref_1.start,
                                                     mod_ref_1.end,
                                                     counter),
                          x))

        unmod, mods = v._current_mods(mod_ref_2.sequence)
        mod_ref_2_seqs = set(remove_nesting(v._apply_mods(unmod, mods)))
        counter = 0
        for x in mod_ref_2_seqs:
            counter += 1
            output_ref.append(
                FASTA_REF('{}:{}-{}_TRANS-{}'.format(mod_ref_2.chromosome,
                                                     mod_ref_2.start,
                                                     mod_ref_2.end,
                                                     counter),
                          x))

        seqs[v.chromosome_1].pop(seqs[v.chromosome_1].index(mod_ref_1))
        seqs[v.chromosome_2].pop(seqs[v.chromosome_2].index(mod_ref_2))
        
        

    for v in vl['CNV']:
        try:
            index_start = search_references(seqs[v.chromosome], v.start_position)
            index_end = search_references(seqs[v.chromosome], v.end_position)
        except IndexError:
            continue

        prev_entries = seqs[v.chromosome][index_start:index_end]

        for pe in prev_entries:
            v.set_relative_start(pe.start)
            counter = 0
            for mod_ref in v.modify(pe.sequence):
                counter += 1
                output_ref.append(
                    FASTA_REF('{}:{}-{}_copy-{}'.format(pe.chromosome,
                                                        pe.start,
                                                        pe.end,
                                                        counter),
                             mod_ref.upper()))

        for pe in prev_entries:
            seqs[v.chromosome].pop(seqs[v.chromosome].index(pe))

    for k in seqs.keys():
        for fasta in seqs[k]:
            output_ref.append(
                FASTA_REF('{}:{}-{}'.format(fasta.chromosome,
                                            fasta.start,
                                            fasta.end),
                          fasta.sequence))

    out_file = open(args.output, 'wb')
    for fasta in output_ref:
        out_file.write('>{}'.format(fasta.name) + '\n')
        sequence = fasta.seq.replace('_', '').replace('X', '')
        
        out_file.write(
            '\n'.join(
                [sequence[X:X+80] for X in range(0, len(sequence), 80)]))
        out_file.write('\n')

    out_file.flush()
    out_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--transvar_variants',
                        help='Formatted variant list exported from TransVar tool')
    parser.add_argument('-r', '--reference',
                        help='Reference gneome',
                        required=True)
    parser.add_argument('-v', '--variants',
                        help="Variants")
    parser.add_argument('-o', '--output',
                        help="Output destination of modified reference",
                        required=True)

    args = parser.parse_args()
    
    main(args)
