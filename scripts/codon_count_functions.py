# codon_count_functions.py
# Functions to calculate codon_usage_table
# As used by DnaOptimizationProblem from DNAchisel
# These functions were inspired by Benjamin Lee's Codon Adaptation Index package:
#   https://github.com/Benjamin-Lee/CodonAdaptationIndex
# Edward Wallace, Edward.Wallace@ed.ac.uk, February 2025

from itertools import chain
from scipy.stats import gmean
from collections import Counter
from Bio import SeqIO
import Bio.Data.CodonTable as ct

# get rid of Biopython warning
import warnings
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)

def _init_codonsbyaa(genetic_code_dict):

    # invert the genetic code dictionary to map each amino acid to its codons
    codons_for_amino_acid = {}
    for codon, amino_acid in genetic_code_dict.items():
        codons_for_amino_acid[amino_acid] = codons_for_amino_acid.get(amino_acid, [])
        codons_for_amino_acid[amino_acid].append(codon)

    # create dictionary of codons by amino acid
    # Example: {'L': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'M': ['ATG']...}
    return codons_for_amino_acid

_codons_by_aa = {
    # dictionary of codons for each amino acid
    k: _init_codonsbyaa(v.forward_table) for k, v in ct.unambiguous_dna_by_id.items()
}


def Count3mers(sequences):
    r"""Counts the 3-mers in frame in a set of sequences.

    Args:
        sequences (list): The  set of sequences.

    Returns:
        dict: Counts of in-frame 3-mers.

    Raises:
        ValueError: When an invalid sequence is provided or a list is not provided.
    """

    if not isinstance(sequences, (list, tuple)):
        raise ValueError(
            "Be sure to pass a list of sequences, not a single sequence. "
        )

    # ensure all input sequences are divisible by three
    for sequence in sequences:
        if len(sequence) % 3 != 0:
            raise ValueError("Input sequence not divisible by three")
        if not sequence:
            raise ValueError("Input sequence cannot be empty")

    # count the number of each codon in the sequences
    sequences = (
        (sequence[i : i + 3].upper() for i in range(0, len(sequence), 3))
        for sequence in sequences
    )
    threemers = chain.from_iterable(
        sequences
    )  # flat list of all codons (to be used for counting)
    counts = Counter(threemers)

    return counts


def CodonFrequencyByAA(sequences, genetic_code=11, outputtype="count"):
    r"""Calculates the synonymous codon frequency by AA for a set of sequences.

    Args:
        sequences (list): The reference set of sequences.
        genetic_code (int, optional): The translation table to use. Defaults to 11, the standard genetic code.
        outputtype: "count" to return codon counts, "freq" to return frequencies proportional to aa.

    Returns:
        dict: For each amino acid, a dict of codon counts or frequencies

    Raises:
        ValueError: When an invalid sequence is provided or a list is not provided.
    """
    
    counts = Count3mers(sequences)
    
    # determine the synonymous codons for the genetic code
    codons_by_aa = _codons_by_aa[genetic_code]

    # hold the result as it is being calculated
    counts_by_aa = {}

    # calculate counts organized by aa
    for aa, codons in codons_by_aa.items():
        miniresult = {}
        for codon in codons:
            # note: check what happens if zero?
            miniresult[codon] = counts[codon]
        counts_by_aa[aa] = miniresult 

    if (outputtype == "count") :
       return counts_by_aa
    elif (outputtype == "freq") :
        freqs_by_aa = {}
        for aa, codons in codons_by_aa.items():
            minifreqs = {}
            miniresult = counts_by_aa[aa]
            total_by_aa = sum(miniresult.values())
            for codon in codons:
                minifreqs[codon] = miniresult[codon] / total_by_aa
            freqs_by_aa[aa] = minifreqs 
        return freqs_by_aa


def fastaFileToStrings(fastafile):
    r"""Helper function to read in fasta file of sequences and convert to plain strings.
    
    Args:
        fastafile (str): A fasta filepath or handle
    
    Returns:
        list: For each record in the fasta file, a string containing the sequence
    """
    
    records = list(SeqIO.parse(fastafile, "fasta"))
    sequencesonly = [str(record.seq) for record in records]
    return sequencesonly
