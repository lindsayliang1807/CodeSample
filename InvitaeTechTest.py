import unittest
import InvitaeTech as main

TEST_INPUT1 = 'input1.tsv'
TEST_INPUT2 = 'input2.tsv'
TEST_OUT = 'out.tsv'
TEST_INPUT1_INVALID = 'input1_invalid.tsv'


class TestQuery(unittest.TestCase):
    '''
    Testing Invitae Bioinformatics Exercise B.3
    '''

    def test_get_transcript_dict_valid(self):
        '''
        Test opening input1 (transcript_file) and storing
        it as a dict
        '''
        expected_dict = {
            'TR1': {
                'gen_chrom': 'CHR1',
                'gen_start': 3,
                'cigar': '8M7D6M2I2M11D7M'},
            'TR2': {
                'gen_chrom': 'CHR2',
                'gen_start': 10,
                'cigar': '20M'}}

        self.assertEqual(
            main.get_transcript_dict(transcript_file=TEST_INPUT1),
            expected_dict)

    def test_get_trans_dict_invalid(self):
        '''
        Test that a transcript file with the incorrect number
        of columns raises an error
        '''
        self.assertRaises(
            ValueError,
            main.get_transcript_dict, TEST_INPUT1_INVALID)

    def test_is_cigar_valid_ok(self):
        '''
        Test is_cigar_valid with valid CIGARs
        '''
        self.assertEqual(
            main.is_cigar_valid(cigar_str='8M7D6M2I2M11D7M'),
            True)
        self.assertEqual(main.is_cigar_valid(cigar_str='20M'), True)

    def test_is_cigar_valid_type(self):
        '''
        Test is_cigar_valid with a CIGAR with invalid qualifiers
        (not M, I, or D)
        '''
        self.assertEqual(main.is_cigar_valid(cigar_str='2M3N7D'), False)

    def test_is_cigar_valid_format(self):
        '''
        Test is_cigar_valid with invalid CIGARs
        '''
        self.assertEqual(main.is_cigar_valid(cigar_str='2M3'), False)
        self.assertEqual(main.is_cigar_valid(cigar_str='M3'), False)
        self.assertEqual(main.is_cigar_valid(cigar_str=''), False)

    def test_cigar_str_to_list(self):
        '''
        Test ciar_str_to_list; converting a cigar string
        to an expanded list of qualifiers (M, D, I)
        '''
        self.assertListEqual(
            main.cigar_str_to_list(cigar_str='1M'),
            ['M'])
        self.assertListEqual(
            main.cigar_str_to_list(cigar_str='2D1M'),
            ['D', 'D', 'M'])
        self.assertListEqual(
            main.cigar_str_to_list(cigar_str='8M7D6M2I2M11D7M'),
            ['M', 'M', 'M', 'M', 'M', 'M', 'M', 'M', 'D',
             'D', 'D', 'D', 'D', 'D', 'D', 'M', 'M', 'M',
             'M', 'M', 'M', 'I', 'I', 'M', 'M', 'D', 'D',
             'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D',
             'M', 'M', 'M', 'M', 'M', 'M', 'M'])

    def test_get_genome_pos_m_offset(self):
        '''
        Test get_genome_pos at an M position with an offset
        (genomic position > 0) before traversing I or D positions
        '''
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=0,
                cigar_str='10M',
                start_pos=3),
            3)
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=4,
                cigar_str='8M7D6M2I2M11D7M',
                start_pos=3),
            7)
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=19,
                cigar_str='8M7D6M2I2M11D7M',
                start_pos=3),
            38)

    def test_get_genome_pos_m_no_offset(self):
        '''
        Test get_genome_pos at an M position without an offset
        (genomic position at 0) before traversing I or D positions
        '''
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=0,
                cigar_str='10M',
                start_pos=0),
            0)
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=4,
                cigar_str='8M7D6M2I2M11D7M',
                start_pos=0),
            4)
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=8,
                cigar_str='8M7D6M2I2M11D7M',
                start_pos=0),
            15)

    def test_get_genome_pos_insert(self):
        '''
        Test get_genome_pos where the query position
        is in an insertion - assumed desired output is an NA
        '''
        # Check position before an insertion
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=4,
                cigar_str='3M2D1I1M1I2M',
                start_pos=2),
            7)
        # Check position after an insertion
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=6,
                cigar_str='3M2D1I1M1I2M',
                start_pos=2),
            8)
        # Check position in the insertion
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=3,
                cigar_str='3M2D1I1M1I2M',
                start_pos=2),
            'NA')
        # Check CIGAR of only insertions
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=1,
                cigar_str='3I',
                start_pos=2),
            'NA')

    def test_get_genome_pos_deletion(self):
        '''
        Test get_genome_pos where the query position
        is in a deletion.
        '''
        # Check the position before the del
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=17,
                cigar_str='8M7D6M2I2M11D7M',
                start_pos=3),
            25)
        # Check the first position after a del
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=8,
                cigar_str='8M7D6M2I2M11D7M',
                start_pos=3),
            18)

    def test_get_genome_pos_combo(self):
        '''
        Test get_genome_pos at an M position with an offset
        after traversing I and D positions
        '''
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=13,
                cigar_str='8M7D6M2I2M11D7M',
                start_pos=3),
            23)
        # Test insertion position next to a deletion
        self.assertEqual(
            main.get_genome_pos(
                transcript_coord=3,
                cigar_str='3M2D1I1M1I2M',
                start_pos=2),
            'NA')

    def test_get_genome_pos_invalid(self):
        '''
        Test get_genome_pos for a query position that's outside
        of the transcript
        '''
        self.assertRaises(
            ValueError,
            main.get_genome_pos, 8, '3M2D1I1M1I2M', 2)

    def test_query_transcript(self):
        '''
        Test query_tramscript method from input files to output files
        '''
        expected_out = [
            ['TR1', '4', 'CHR1', '7'],
            ['TR2', '0', 'CHR2', '10'],
            ['TR1', '13', 'CHR1', '23'],
            ['TR2', '10', 'CHR2', '20']
        ]

        main.query_transcript(
            transcript_file=TEST_INPUT1,
            query_file=TEST_INPUT2,
            out=TEST_OUT)

        out_list = list()
        with open(TEST_OUT, 'r') as result_f:
            for line in result_f:
                out_list.append([i.strip() for i in line.split('\t')])
        self.assertEqual(out_list, expected_out)

if __name__ == '__main__':
    unittest.main()
