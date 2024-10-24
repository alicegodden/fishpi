import unittest
from FishPi import classify_te_type, introduce_mismatches

class TestFishPi(unittest.TestCase):

    def test_classify_te_type(self):
        self.assertEqual(classify_te_type("hat"), "DNA")
        self.assertEqual(classify_te_type("gypsy"), "LTR")
        self.assertEqual(classify_te_type("unknown"), "UNKNOWN")

    def test_introduce_mismatches(self):
        seed = "ATGCATGC"
        result = introduce_mismatches(seed, 2)
        self.assertEqual(len(result), len(seed))
        self.assertNotEqual(result, seed)  # With mismatches, it should not be the same

if __name__ == '__main__':
    unittest.main()
