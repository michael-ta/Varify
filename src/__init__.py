import abc
import random
import re
import collections
import string

# Define a new variant to be added into the synthetic reference
# This is the base class for defining variants. Define the class
# variable ATTRIBUTES in implementation. Each attribute mapping
# should be a tuple similar to (name, type_initializer)
# Input values will be set to the type during initialization


# Define a named tuple for modifications
MOD = collections.namedtuple('MOD', ['type', 'index', 'effect'])

# strip unused / blank fields before initializing
class VariantDefinition(object):
    ATTRIBUTES = []
    TRANSLATION_TABLE = {}
    TRANS_COUNT = 0
    def __init__(self, *args, **kwargs):
        self.type = kwargs['var_type']
        self._rel_start = None
        self._map_arguments_to_attributes(*args)

    @abc.abstractmethod
    def modify(self, reference):
        pass

    def set_relative_start(self, position):
        self._rel_start = position

    def _map_arguments_to_attributes(self, *args):
        pairs = zip(self.ATTRIBUTES, args)
        for k,v in pairs:
            setattr(self, k[0], k[1](v))

    def _apply_mods(self, reference, mods, tags_only=False):
        modified_reference = ''
        prev_index = 0
        sorted_mods = sorted(mods, key=lambda x:x.index)
        for i, m in enumerate(sorted_mods):
            if m.type == 'INS':
                modified_reference += reference[prev_index:m.index] + '_{}_'.format(m.effect)
                prev_index = m.index
            elif m.type == 'TRANS':
                # this returns it, but we really want to apply the modification
                modified_reference += reference[prev_index:m.index]
                if not tags_only:
                    modified_reference += self._get_trans_reference(VariantDefinition.TRANSLATION_TABLE[m.effect], m.effect)

                    return (self._apply_mods(*self._current_mods(modified_reference)),
                    self._apply_mods('X' * m.index + reference[m.index:], sorted_mods[i+1:]))
                else:
                    # not tag only so append the MOD to the reference
                    modified_reference += '|({})|'.format(m.effect)
                
                prev_index = m.index
            elif m.type =='INV':
                modified_reference += reference[prev_index:m.index] 
                if not tags_only:
                    inv_seq, other_seq = self._get_inv_reference(
                        VariantDefinition.TRANSLATION_TABLE[m.effect], m.effect)

                    modified_reference += inv_seq

                   
                    return (self._apply_mods(*self._current_mods(modified_reference)),
                    self._apply_mods(*self._current_mods('X' * m.index + reverse_complement(reference[m.index:]) + other_seq)))
                    #self._apply_mods('X' * m.index + reference[m.index:], sorted_mods[i+1:]))

                else:
                    modified_reference += '|({})|'.format(m.effect)
                prev_index = m.index

        return modified_reference + reference[prev_index:]
            
    def _get_trans_reference(self, reference, effect):
        count, event = effect.split('-')
        if int(count) % 2 == 0:
            query = '|({}-{})|'.format(int(count) - 1, event)
        else:
            query = '|({}-{})|'.format(int(count) + 1, event)

        index = reference.find(query)
        return reference[index + len(query):]

    def _get_inv_reference(self, reference, effect):
        count, event = effect.split('-')
        if int(count) % 2 == 0:
            query = '|({}-{})|'.format(int(count) - 1, event)
        else:
            query = '|({}-{})|'.format(int(count) + 1, event)
        index = reference.find(query)
        return reverse_complement(reference[:index]), reference[index + len(query):]
        
    def _current_mods(self, reference):
        ''' if other insertions have been made we will need to strip them from the reference
            before we proceed
        '''
        # split out the insertions
        data = reference.split('_')
        mods = []
        index = 0
        unmodified_ref = ''
        for x in range(len(data)):
            if x % 2 == 0:
                # this is part of the unmodified reference
                data_trans = data[x].split('|')
                if len(data_trans) > 1:
                    for trans_index in range(len(data_trans)):
                        if '(' not in data_trans[trans_index]:
                            index += len(data_trans[trans_index])
                            unmodified_ref += data_trans[trans_index]
                        else:
                            effect = data_trans[trans_index][1:-1]
                            mods.append(
                                MOD(effect.split('-')[1],
                                    index,
                                    effect))
                else:
                    index += len(data[x])
                    unmodified_ref += data_trans[0]
            else:
                index += len(data[x])
                
                mods.append(
                    MOD('INS',
                        index,
                        data[x]))

        return unmodified_ref, mods


    def _identify_translocations(self, ref):
         positions = self._identify_index(ref)
         events = re.findall('\|\(([^\|]*)\)', ref)

         return zip(positions, events)

    def _indentify_index(self, ref):
        positions = []
        for i,x in enumerate(ref):
            if x == '|':
                positions.append((i,x))
        return positions


class SubstitutionVariantDefinition(VariantDefinition):
    ''' This class defines a single base variant
    '''
    ATTRIBUTES = [('chromosome', str),
                  ('position', int),
                  ('ref', str),
                  ('alt', str)]

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)

    def modify(self, reference):
        # no modifications should be made to this base as of yet...
        # convert from 1 - base indexing to 0 base
        index = self.position - self._rel_start - 1 
        assert(reference[index] == self.ref)
        start_seq = ''
        end_seq = ''
        try:
            start_seq = reference[:index]
        except IndexError:
            raise Exception("Index for variant does not exists in reference")
        try:
            end_seq = reference[index + 1:]
        except IndexError:
            # if it didn't error out from before, it means that the position of
            # the variant appears at the very end of the read
            end_seq = ''
        return start_seq + self.alt + end_seq


class DeletionVariantDefinition(VariantDefinition):
    ATTRIBUTES = [('chromosome', str),
                  ('start_position', int),
                  ('end_position', int)]

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        # the actual position / index of the deleted base should be + 1
        self.start_position += 1

    def modify(self, reference):
        # a couple of cases should be considered here for the deletion
        # 1. complete overlap
        #   when the deletion appears before the start of reference and ends after
        # 2. complete overlap of one end
        #   when the deletion appears before the reference and ends in the middle
        #   or the deletion appears in the middle of the refernces and ends
        #   after
        # 3. contained within the reference
        #   when the deletion appears somewhere in the middle of the reference
        #   and ends before the end
        len_reference = len(reference)
        modified_reference = ''

        if self.start_position <= self._rel_start:
            # subtract 1 from the index possibly
            if self.end_position >= self._rel_start + len_reference - 1:
                # handle case 1
                # don't do anything, an empty string should be returned
                deletion_start_index = 0
                deletion_end_index = len_reference
            else:
                # handle case 2
                # deletion begins at or before the start of the sequence and
                # terminates in the middle
                deletion_start_index = 0
                deletion_end_index = self.end_position - self._rel_start + 1

        elif self.end_position >= self._rel_start + len_reference - 1:
            # handle case 2
            # deletion begins in the middle of the sequence and terminates at the end
            deletion_start_index = self.start_position - self._rel_start
            deletion_end_index = len_reference
            pass
        else:
            # handle case 3
            deletion_start_index = self.start_position - self._rel_start
            deletion_end_index = self.end_position - self._rel_start + 1
        modified_reference += reference[:deletion_start_index]
        modified_reference += 'X' * (deletion_end_index - deletion_start_index)
        # why does this have to be incremented by 1
        modified_reference += reference[deletion_end_index:]

        return modified_reference 


class InsertionVariantDefinition(VariantDefinition):
    ''' Insertions can have predefined bases or simply define the number of bases
        that are inserted, at which point the inserted bases will be randomly
        generated
    '''
    ATTRIBUTES = [('chromosome', str),
                  ('position', int),
                  ('length', int),
                  ('bases', str)]

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        try:
            self.bases
        except AttributeError:
            self._randomize_bases()

    def modify(self, reference):
        # a couple of cases should be considered here for insertions
        # 1. extension before
        #   when the insertion occurs at the base right before the start of
        #   reference
        # 2. extension after
        #   when the insertion occurs at the base right after the end of the
        #   reference
        # 3. contained within the reference
        #   when the insertion appears somewhere in the middle of the reference
        # 4. preexisting insertion
        #   when there is a preexisting insertion to worry about; we have to split
        #   the reference read and the piece it back together with all inserted bases

        # no modifications necessary
        if self.position < self._rel_start - 1:
            return reference
        elif self.position > self._rel_start + len(reference):
            return reference

        reference, mods = self._current_mods(reference)
        
        len_reference = len(reference)

        # handle case 1
        if self.position == self._rel_start - 1:
            insertion_start_index = 0
        # handle case 2
        elif self.position == self._rel_start + len(reference):
            insertion_start_index = len_reference
        # handle case 3
        else:
            insertion_start_index = self.position - self._rel_start + 1

        mods.append(MOD(
                    self.type,
                    insertion_start_index,
                    self.bases))

        return self._apply_mods(reference, mods)

    '''
    def _apply_mods(self, reference, mods):
        modified_reference = ''
        prev_index = 0
        for index, bases in sorted(mods, key=lambda x:x[0]):
            modified_reference += reference[prev_index:index] + '_{}_'.format(bases)
            prev_index = index

        return modified_reference + reference[prev_index:]
    
     
    def _current_mods(self, reference):
        #if other insertions have been made we will need to strip them from the reference
        #    before we proceed

        data = reference.split('_')
        mods = []
        unmodified_ref = ''.join([data[x] for x in range(0, len(data), 2)])
        pos = 0
        for x in range(1, len(data), 2):
            pos += len(data[x-1])
            mods.append((pos, data[x]))

        return unmodified_ref, mods
    '''

    def _randomize_bases(self):
        ''' if a set of bases where not passed into the constructor for the insertion,
            then generate set of randomized bases for use
        '''
        
        self.bases = ''.join([random.choice('AGCT') for x in range(self.length)])


class TranslocationVariantDefinition(VariantDefinition):
    ATTRIBUTES = [('chromosome_1', str),
                  ('position_1', int),
                  ('chromosome_2', str),
                  ('position_2', int)]

    # maintain a dictionary containing all translocation reads
    # format will follow as defined
    # {EVENT-ID: SEQUENCE, ... }


    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.ids = []
        self._set_event_ids()

    def modify(self, reference1, reference2):
        ''' given two reference sequences, return the translocated sequences
        '''
        modified_reference1 = ''
        modified_reference2 = ''

        if self.type == 'INV':
            index1 = self.position_1 - self._rel_start_1 + 1
            index2 = self.position_2 - self._rel_start_2 + 1

            unmod_ref_1, mods_1 = self._current_mods(reference1)
            unmod_ref_2, mods_2 = self._current_mods(reference2)

            mods_1.append(MOD(self.type,
                              index1,
                              self.ids[0]))
            mods_2.append(MOD(self.type,
                              index2,
                              self.ids[1]))

            modified_reference1 = self._apply_mods(unmod_ref_1, mods_1, tags_only=True)
            modified_reference2 = self._apply_mods(unmod_ref_2, mods_2, tags_only=True)

        elif self.type == 'TRANS':
            index1 = self.position_1 - self._rel_start_1 + 1
            index2 = self.position_2 - self._rel_start_2 + 1
            unmod_ref_1, mods_1 = self._current_mods(reference1)
            unmod_ref_2, mods_2 = self._current_mods(reference2)

            mods_1.append(MOD(self.type,
                               index1,
                               self.ids[0]))
            mods_2.append(MOD(self.type,
                               index2,
                               self.ids[1]))

            modified_reference1 = self._apply_mods(unmod_ref_1, mods_1, tags_only=True)
            modified_reference2 = self._apply_mods(unmod_ref_2, mods_2, tags_only=True)

        
        self._update_trans_table(modified_reference1)
        self._update_trans_table(modified_reference2)

        return modified_reference1, modified_reference2

    def _insert_delimiter(self, reference, index, event_id):
        ''' given the reference and the index position, insert into the reference
            the '|' delimiter for structural variants
        '''
        data = reference.split('_')
        current_index = 0
        for x in range(len(data)):
            if current_index + len(data[x]) > index:
                data[x] = data[x][:index - current_index] + '|({})|'.format(event_id) + data[x][index-current_index:]
                break
            if x % 2 == 0:
                current_index += len(data[x])
        return '_'.join(data)

    def _set_event_ids(self):
        ''' set the id for this translocation event
        '''
        VariantDefinition.TRANS_COUNT += 1
        self.ids.append('{}-{}'.format(VariantDefinition.TRANS_COUNT, self.type))
        VariantDefinition.TRANS_COUNT += 1
        self.ids.append('{}-{}'.format(VariantDefinition.TRANS_COUNT, self.type))

    def _update_trans_table(self, ref):
        ''' updated a reference or add a new reference to the translation table
            remember that the event_id needs to map to the bases that are going to be
            used to modify the reference -- therefore even ids are - 1 and odd are + 1
        '''
        for mod in self._current_mods(ref)[1]:
            if mod.type in ['TRANS', 'INV', 'RTRANS', 'RINV']:
                count, event = mod.effect.split('-')
                if int(count) % 2 == 0:
                    VariantDefinition.TRANSLATION_TABLE['{}-{}'.format(int(count) - 1, event)] = ref
                else:
                    VariantDefinition.TRANSLATION_TABLE['{}-{}'.format(int(count) + 1, event)] = ref
        
    def merge_mods(self, reference1, reference2):
        ''' merge together 2 different translocation modifications
        '''

    def set_relative_start(self, position1, position2):
        ''' set the relative start positions of the references
            position 1 should correspond to the refernece1
            position 2 should correspond to the reference2
        '''
        self._rel_start_1 = position1
        self._rel_start_2 = position2


class CopynumberVariantDefinition(VariantDefinition):
    ''' can only define copy number gains at the moment
    '''
    ATTRIBUTES = [('chromosome', str),
                  ('start_position', int),
                  ('end_position', int),
                  ('copy_number', int)]

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)

    def modify(self, reference):
        len_reference = len(reference)
        
        if self.start_position <= self._rel_start:
            # subtract 1 from the index possibly
            if self.end_position >= self._rel_start + len_reference - 1:
                # handle case 1
                # the CNV completely encompsses the reference
                cnv_start_index = 0
                cnv_end_index = len_reference
            else:
                # handle case 2
                # CNV begins at or before the start of the sequence and
                # terminates in the middle
                cnv_start_index = 0
                cnv_end_index = self.end_position - self._rel_start + 1

        elif self.end_position >= self._rel_start + len_reference - 1:
            # handle case 2
            # cnv begins in the middle of the sequence and terminates at the end
            cnv_start_index = self.start_position - self._rel_start
            cnv_end_index = len_reference
            pass
        else:
            # handle case 3
            cnv_start_index = self.start_position - self._rel_start
            cnv_end_index = self.end_position - self._rel_start + 1
        copy_reference = [reference[cnv_start_index:cnv_end_index] for x in range(self.copy_number - 1)]

        return [reference] + copy_reference 


class Reference(object):
    def __init__(self):
        # modification list
        self.mod_list = []

    def get_reference(self):
        ''' for each item in the mod_list, apply the modification to the
            reference sequence accordingly and return the modified reference
            sequence
            TODO: cache reference for quicker access?
        '''
        pass


def reverse_complement(read):
    trans_table = {'A': 'T',
                   'T': 'A',
                   'C': 'G',
                   'G': 'C'}
    complement = ''
    event_trans = False
    event_ins = False
    event_string = ''
    
    for b in range(1, len(read) + 1):
        if read[-b] == '|':
            if event_trans:
                event_trans = False
                complement += event_string
                event_string = ''
            else:
                event_trans = True
            complement += read[-b]
            continue
        elif read[-b] == '_':
            if event_ins:
                event_ins = False
            else:
                event_ins = True
            complement += read[-b]
            continue
        else:
            if event_trans:
                event_string = read[-b] + event_string
            else:
                try:
                    complement += trans_table[read[-b]]
                except KeyError:
                    complement += read[-b]
    return complement

def remove_nesting(v):
    ''' remove nested tuples
    '''
    if type(v) == tuple:
        try:
            v[1]
            return [v[0]] + remove_nesting(v[1])
        except IndexError:
            return [v[0]]
    else:
        return [v]

VARIANT_CLASS_TABLE = {
    'SNP': SubstitutionVariantDefinition,
    'TRANS': TranslocationVariantDefinition,
    'INV': TranslocationVariantDefinition,
    'INS': InsertionVariantDefinition,
    'DEL': DeletionVariantDefinition,
    'CNV': CopynumberVariantDefinition}
