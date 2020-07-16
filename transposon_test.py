#Before running transposon_test, in transposon_finder line 72 set transposon_sim = 1
import unittest
from transposon_finder import window_compiler, window_streamliner, transposon_identifier, seq_2_order
test_seq_1 = "AAAGGGCCCTTT"
test_seq_2 = "CCCAAATTTGGG"



class finder_test(unittest.TestCase):
	def test_window_library(self):
		run = window_streamliner(test_seq_1, test_seq_2)
		expect = [[0, 3], [3, 6], [6, 9], [9, 12]]
		self.assertEqual(run, expect)


def test_transposon_identifier(self):
		run = transposon_identifier(test_seq_1,test_seq_2)
		expect = ["CCC", "AAA", "TTT", "GGG"]
		self.assertEqual(run, expect)

		
def seq_2_order(test_seq_1, test_seq_2):
		run = seq_2_order(test_seq_1, test_seq_2)
		expect = {'CCC': 1, 'AAA': 2, 'TTT': 3, 'GGG': 4}
		self.assertEqual(run, expect)

	
def seq_1_order(test_seq_1, test_seq_2):
		run = seq_1_order
		expect = [2, 4, 1, 3]
		self.assertEqual(run, expect)



unittest.main()