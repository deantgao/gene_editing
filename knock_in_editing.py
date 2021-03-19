#!/usr/bin/env python
# coding: utf-8

# This is a straightforward implementation to design HA primers from ENST transcript ID. It does not account for potential nuances in how Cas9, gene knock ins, and homologous arm binding might behave. Moreover, it identifies ALL candidates templates that meet the below specified requirements, but it does not yet weight each candidate based on its likelihood of successful integration. This program would need, in the very least, a better weighting algorithm before deployment.
# 
# It is important to acknowledge the assumptions that went into this design, which may or may not meet all functional requirements.
# 
# These assumptions include:
# - PAM (NGG) sequence is used to identify candidate cut/insertion sites.
# - ~50 nt for each homologous arm, possibly extended on the 3' arm to ensure G/C clamp.
# - Left homologous arm directly precedes the beginning of PAM sequence; right homologous arm is offset by 5 nt AFTER the beginning of PAM sequence.
# - Left homologous arm is forward encoded; right homologous arm is reverse complement encoded.
# - Exactly 3 G or C nucleotides among the last 5 nucleotides of the sequence template.
# - No more than 4 consecutive runs of a single base per arm.
# - No more than 4 consecutive di-nucleotide repeats per arm.
# - 40-60% GC content in each arm.
# - Melting temperature between 55-80 celsius in each arm.
# - No more than a 3 degree celsius difference between left and right arms.

# In[1]:


import os
import sys
import requests
import itertools
import ensembl_rest
from pprint import pprint


# In[2]:


SERVER = "http://rest.ensembl.org" 
GFP = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCTGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGCGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA'


# In[3]:


def genQuery(transcript_id):
    '''
    Returns --> the sequence associated with a unique Ensembl transcript ID
    '''
    seqext = f'/sequence/id/{transcript_id}?'
    return seqext
    
def fetchEndpoint(server, query, content_type):
    '''
    Retrieves the specified API endpoint from Ensembl.
    '''
    result = requests.get(server + query, headers = { "Content-Type" : f'{content_type}'})
    if not result:
        result.raise_for_status()
        sys.exit()
    if content_type == 'application/json':
        return result.json()
    return result.text

def findCandidatePAMs(sequence, arm_length):
    '''
    Parses the nucleotide sequence and identifies the start index of any candidate PAM ('NGG').
    
    Returns --> in order list of the beginning indices of all candidate PAM
    '''
    start_idx = []
    for i in range(len(sequence) - 1):
        if i > arm_length and sequence[i] == sequence[i + 1] == 'G':
            start_idx.append(i - 1)
    return start_idx

def findBeginSeqIndices(sequence, arm_length):
    '''
    Returns --> list of indices, each one indicating the starting index of the 5' homologous arm
    '''
    result = findCandidatePAMs(sequence, arm_length)
    indices = []
    for idx in result:
        indices.append(idx - arm_length)
    return indices

def findHomologousArms(sequence, idx, arm_length):
    '''
    Given the start idx and arm length, finds the nucleotide sequences to generate the forward and reverse primers
    The forward primer is the 50 nt sequence preceding PAM. The reverse primer is the ~50 nt sequence following
    PAM, but NOT including PAM.
    
    Returns --> fwd primer, reverse primer
    '''
    def handleEndingGC(sequence, primer, idx):
        '''
        Make sure the reverse 3' primer ends in G or C.
        '''
        ptr = idx
        while sequence[ptr] != 'G' and sequence[ptr] != 'C':
            ptr -= 1
            primer = sequence[ptr] + primer
        return primer
    
    nucleo_mapping = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    start = idx
    PAM_offset = 5
    fwd_primer = sequence[start : start + arm_length]
    temp_primer = sequence[start + arm_length + PAM_offset: start + (2 * arm_length) + PAM_offset]
    rev_primer = ''
    for nucleo in temp_primer:
        rev_primer += nucleo_mapping[nucleo]
    rev_primer = handleEndingGC(sequence, rev_primer, start + arm_length + PAM_offset)
    rev_primer = rev_primer[::-1]
    return (fwd_primer, rev_primer)

def calcNucleoContent(sequence):
    '''
    Computes the count of G/C and A/T in the sequence.
    '''
    num_gc = sequence.count('G') + sequence.count('C')
    num_at = len(sequence) - num_gc
    return (num_gc, num_at)

def calcMeltingTm(GC_content, AT_content):
    '''
    Calculates the melting temperatire using formula:
    Tm = 64.9 + 41 * (yG + zC - 16.4) / (wA + xT + yG + zC)
    
    Returns --> Computed melting temperature
    '''
    melt_tm = 64.9 + 41 * (GC_content - 16.4) / (AT_content + GC_content)
    return melt_tm

def diffInTmBtwnArms(tm_arm1, tm_arm2):
    '''
    The maximum difference should be no more than 3 celsius between forward and reverse primer.
    
    Returns --> Difference in temperature between fwd and rev primer
    '''
    return abs(tm_arm1 - tm_arm2) 

def confirmGCon3Prime(rev_primer):
    '''
    Counts number of G/C within the last 5 nt of the sequence.
    
    Returns --> Count of G/C 
    '''
    gc_on3 = rev_primer[-5 : ].count('G') + rev_primer[-5 : ].count('C')
    return gc_on3

def isSeqValid(primer):
    '''
    Determines whether a single arm primer is valid
    
    Returns --> True if primer is valid else False
    '''
    def moreThan4Runs(sequence):
        '''
        Checks whether there is any single base that appears more than 4 times consecutively in the sequence.

        Returns --> True if there are more than 4 else False
        '''
        prev = sequence[0]
        consecutive_runs = 1
        for nucleo in sequence[1 : ]:
            if prev == nucleo:
                consecutive_runs += 1
            else:
                prev = nucleo
                consecutive_runs = 1
            if consecutive_runs > 4:
                return True
        return False
    
    def moreThan4Repeats(sequence):
        '''
        Checks whether there are any di-nucleotide sequences that repeat more than 4 times consecutively in 
        the sequence.
        
        Returns --> True if there are more than 4 else False
        '''
        nucleos = 'ATCG'
        repeats = set()
        for pair in itertools.combinations(nucleos, 2):
            repeats.add(''.join(pair * 5))
        for repeat in repeats:
            if repeat in sequence:
                return True
        return False
    
    if moreThan4Runs(primer):
        return False
    if moreThan4Repeats(primer):
        return False
    return True

def isPrimerPairValid(cut_idx, fwd, rev):
    '''
    Determines whether two homologous arms combined are valid.
    
    Returns --> True if pair is valid else False
    '''
    fwd_gc_content, fwd_at_content = calcNucleoContent(fwd)
    fwd_percent_gc = fwd_gc_content / len(fwd)
    rev_gc_content, rev_at_content = calcNucleoContent(rev)
    rev_percent_gc = rev_gc_content / len(rev)
    fwd_tm = calcMeltingTm(fwd_gc_content, fwd_at_content)
    rev_tm = calcMeltingTm(rev_gc_content, rev_at_content)
    temp_diff = diffInTmBtwnArms(fwd_tm, rev_tm)
    if confirmGCon3Prime(rev) != 3:
        return False
    if fwd_percent_gc < 0.4 or fwd_percent_gc > 0.6 or rev_percent_gc < 0.4 or rev_percent_gc > 0.6:
        return False
    if not isSeqValid(fwd) or not isSeqValid(rev):
        return False
    if fwd_tm < 55 or fwd_tm > 80 or rev_tm < 55 or rev_tm > 80:
        return False
    if temp_diff > 3:
        return False
    print(f'This is a candidate for exogenous template.')
    print(f'The index for cut site: {cut_idx}')
    print(f'The temp in left arm: ~{int(fwd_tm)} and right arm: ~{int(rev_tm)}')
    print(f'The GC content of left arm: {round(fwd_gc_content / len(fwd), 2)} and right arm: {round(rev_gc_content / len(rev), 2)}')
    print(f'The template ends in {rev[-1]} on the 3\' end and has a total 3 G/C nucleotides in the last 5 nts of the sequence ({rev[-5 : ]})\n')
    return True

def genTemplate(fwd_primer, rev_primer, GFP):
    '''
    Generates candidate template by inserting the GFP sequence between the left and right homologous arms.
    
    Returns --> string of exogenous template consisting of the GFP sequence (in lowercase for easy 
    distinguishing between GFP and homologous arms) flanked by fwd homologous arm on left and rev homologous
    arm on right
    ''' 
    return fwd_primer + GFP.lower() + rev_primer


# In[4]:


def genAllCandidateTemplates(transcript_ID, arm_length):
    '''
    Generates all valid exogenous templates flanked by homologous arms.
    Prints a robust profile for each valid template.
    
    Candidate template meets following requirements:
    A sequence whose left homologous arm is the 50 nt preceding PAM and whose right homologous arm is the ~50 
    nt sequence beginning 5 nt AFTER beginning of PAM sequence.
    A forward and reverse complementary primer to flank either end of the GFP sequence.
    No less than 3 and no more than 3 G/C in the last 5 nt on the 3' arm.
    No more than 4 consecutive runs of a single base per arm.
    40-60% GC content in each primer arm.
    Melting temperature between 55-80 celsius in each arm.
    No more than a 3 degree celsius difference between left and right arms.
    
    Returns --> (index of cut site (right before PAM), exogenous template sequence with homologous arms)
    '''
    query = genQuery(transcript_ID)
    target_sequence = fetchEndpoint(SERVER, query, 'text/plain')
    candidate_idx = findBeginSeqIndices(target_sequence, arm_length)
    candidate_templates = []
    for idx in candidate_idx:
        fwd, rev = findHomologousArms(target_sequence, idx, arm_length)
        if not isPrimerPairValid(idx + arm_length, fwd, rev):
            continue
        template = genTemplate(fwd, rev, GFP)
        candidate_templates.append((idx + arm_length, template))
    print(f'There are a total {len(candidate_templates)} candidate templates for GFP insertion within the sequence associated with {transcript_ID}')
    return candidate_templates


# In[5]:


transcript_id = 'ENST00000369985'
templates = genAllCandidateTemplates(transcript_id, 50)


# In[6]:


print(f'The candidate templates are as follows:\n{templates}')

