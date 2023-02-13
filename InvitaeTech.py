'''
Invitae Bioinformatics Exercise B.3 (transcript query)
Author: Lindsay Liang

Assumptions:
    - Valid chromosome/position pairs in input1
    - if a query position is at an insertion position in the
      transcript, return NA
    - Only handles M, I, D CIGAR operators
      (not N, for RNA CIGARs or clipping)
    - Empty CIGAR strings are invalid
    - CIGAR with only deletions are invalid

Strengths:
    - For each query in input2, runtime is O(n) where
      n is the length of the transcript (not the cigar), and space efficiency
      is also linear on transcript length
    - Checks user input for correct number of columns, validates CIGAR string

Weaknesses:
    - Processes each query in sequence; could be done in parallel
      after reading in all of input2 and storing it in a dict
    - Could be O(n) where n is the length of the cigar, if using the following
      in get_genome_pos:
        ```
        transcript_idx = -1
        genome_pos = start_pos -1
        for idx in range(0, len(cigar_str)):
            if transcript_idx >= transcript_coord:
                genome_pos = genome_pos + transcript_coord - transcript_idx
                break
            char = cigar_str[idx]
            if char.isnumeric():
                num_bases += char
            elif char == 'M':
                transcript_idx += int(num_bases)
                genome_pos += int(num_bases)
                num_bases = str()
            elif char == 'D':
                genome_pos += int(num_bases)
                num_bases = str()
            elif char == 'I':
                # insertion
                transcript_idx += int(num_bases)
                num_bases = str()
        ```
      But this doesn't handle the case where the query position is in an
      insertion in the transcript. So depending on the desired behavior of
      handling that case, could be more efficient.  This solution would also
      be more space efficient since the expanded cigar string would not have to
      be stored.  However I think that the presented solution in get_genome_pos
      is more intuitive/readable, with the same time complexity.

Testing:
    - InvitaeTechTest.py for unittests.
'''

from argparse import ArgumentParser
import csv
import re


def get_transcript_dict(transcript_file):
    '''
    Builds a dict from the transcript_file key'd by transcript IDs.

    Args:
        transcript_file: path to the tsv of file transcripts.
            4 columns: transcript ID, genomic chromosome, genomic position, and
            cigar string.
    Returns:
        A dict of dicts of transcripts, their chromosome location,
        start position and cigar string.
    Raises:
        ValueError: if there are duplicate transcript ID's.
    '''

    transcripts = dict()
    with open(transcript_file, 'r') as tr_f:
        for line in tr_f:
            fields = [i.strip() for i in line.split('\t')]
            if len(fields) != 4:
                raise ValueError(f'Incorrect format; only {len(fields)} cols.')
            tr_id, chrom, pos, cigar = fields
            if not is_cigar_valid(cigar):
                raise ValueError(f'Invalid cigar {cigar}.')
            if tr_id not in transcripts:
                transcripts[tr_id] = {
                    'gen_chrom': chrom,
                    'gen_start': int(pos),
                    'cigar': cigar
                }
            else:
                # Assume having duplicate id's is user error
                raise ValueError(f'Transcript {tr_id} is not unique.')
    return transcripts


def is_cigar_valid(cigar_str):
    '''
    Validates cigar_str
    Args:
        cigar_str: cigar string to validate
    Returns:
        True or False, depending on if cigar_str is valid
    Raises:
        None
    '''

    return bool(re.match(r'^([0-9]+[M|I|D])+$', cigar_str))


def cigar_str_to_list(cigar_str):
    '''
    Expands the CIGAR string to a list where each character
    is either M|D|I matching the position described in cigar_str.
    Eg 3M1D2I -> [M, M, M, D, I, I]

    Args:
        cigar_str: cigar string
    Returns:
        A list of characters for every position in the CIGAR

    Raises:
        None
    '''
    cigar_list = list()
    num_str = str()
    for char in cigar_str:
        if char.isnumeric():
            # append characters that are numbers for multiple digits
            num_str += char
        else:
            # char.isalpha():
            cigar_list.extend(int(num_str) * [char])
            # reset num_str
            num_str = str()

    return cigar_list


def get_genome_pos(transcript_coord, cigar_str, start_pos):
    '''
    Calculates the genomic position of a transcript coordinate.

    Args:
        transcript_coord: the query coordinate of the transcript/cigar_str
        cigar_str: the properly formated cigar string of the transcript
        start_pos: the genomic start position of the transcript
    Returns:
        genome_pos: the genomic position of the transcript coordinate, or NA if
            the transcript coordinate is in an insertion of the genome
    Raises:
        ValueError: if transcript coord is outside of the length of the cigar
    '''

    cigar_list = cigar_str_to_list(cigar_str)
    genome_pos = start_pos - 1
    transcript_idx = - 1

    # Check if query coordinate is within transcript
    transcript_len = len(cigar_list) - cigar_list.count('D')
    if transcript_coord >= transcript_len:
        raise ValueError(
            f'Query coord {transcript_coord} is outside of the query.')

    for pos in cigar_list:
        if transcript_idx < transcript_coord:
            if pos == 'M':
                # match, both pointers move forward
                transcript_idx += 1
                genome_pos += 1
            elif pos == 'D':
                # tr deletion, only genome_pos pointer moves forward
                genome_pos += 1
            elif pos == 'I':
                # tr insertion, only tr_coord pointer moves forward
                transcript_idx += 1
                # if the transcript coordinate is in an insertion, return NA
                if transcript_idx == transcript_coord:
                    return 'NA'
        else:
            break

    return genome_pos


def query_transcript(transcript_file, query_file, out):
    '''
    For every query in query_file, output a line in out where each line is
    4 cols of the query transcript, the query transcript position, the genomic
    chromosome, and the genomic position of the query position.

    Args:
        transcript_file: path to the tsv of file transcripts.
            4 columns: transcript ID, genomic chromosome, genomic position, and
            cigar string.
        query_file: path to the tsv of queries.
            2 colunns: transcript ID to query, and transcript coordinate
        out: path to output tsv.
            A line for every line in query_file, and 4 coluns: transcript ID,
            transcript coordinate, genomic chrom, and genomic position.
    Returns:
        None; writes an output file.
    Raises:
        ValueError: if either input tsvs have the wrong number of columns
    '''

    # Build transcript dict
    transcripts = get_transcript_dict(transcript_file=transcript_file)

    # Iterate through queries
    with open(out, 'w') as out_f:
        out_writer = csv.writer(out_f, delimiter='\t')
        with open(query_file, 'r') as query_f:
            for line in query_f:
                fields = [i.strip() for i in line.split('\t')]
                if len(fields) != 2:
                    raise ValueError(f'Incorrect format; {len(fields)} cols.')
                tr_id, tr_coord = fields[0], int(fields[1])

                # Get query transcript info from transcript dict
                gen_chrom = transcripts[tr_id]['gen_chrom']
                gen_start = transcripts[tr_id]['gen_start']
                cigar_str = transcripts[tr_id]['cigar']
                gen_pos = get_genome_pos(
                    transcript_coord=tr_coord,
                    cigar_str=cigar_str,
                    start_pos=gen_start)

                out_writer.writerow([tr_id, tr_coord, gen_chrom, gen_pos])
    return


def parse_args():
    '''
    Parse arguments at runtime.
    Args:
        None; takes user input.
    Returns:
        Arguments inputted by user for transcript_fie, query_file, and out.
    Raises:
        None
    '''
    parser = ArgumentParser()
    parser.add_argument('-t', '--transcript-file', type=str, required=True)
    parser.add_argument('-q', '--query-file', type=str, required=True)
    parser.add_argument('-o', '--out', type=str, required=True)

    args = parser.parse_args()
    return args.transcript_file, args.query_file, args.out


def main(transcript_file, query_file, out):
    '''
    Main function, wraps query_transcript

    Args:
        transcript_file: path to the tsv of file transcripts.
            4 columns: transcript ID, genomic chromosome, genomic position, and
            cigar string.
        query_file: path to the tsv of queries.
            2 colunns: transcript ID to query, and transcript coordinate
        out: path to output tsv.
            A line for every line in query_file, and 4 coluns: transcript ID,
            transcript coordinate, genomic chrom, and genomic position.
    Returns:
        None; writes an output file.
    '''
    query_transcript(
        transcript_file=transcript_file,
        query_file=query_file,
        out=out)

if __name__ == '__main__':
    main(
        *parse_args()
    )
