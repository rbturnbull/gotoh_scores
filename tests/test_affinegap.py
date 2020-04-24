# -*- coding: utf-8 -*-
import unittest
import gotoh_scores
import Levenshtein



    

class GotohScoresTest(unittest.TestCase):
    def setUp(self):
        self.scores = gotoh_scores.scores
    
    def levenshtein( self, string1, string2 ):
        g = self.scores( string1, string2 )
        d = Levenshtein.distance( string1, string2 )
    
        print("Testing '%s' and '%s' compared with Levenshtein Distance [%d == %d]" % (string1, string2, -g[0], d) )
        self.assertAlmostEqual( -g[0], d )    
    
    def assertAlmostEqualGotoh(self, g, true_values, delta=0.001):
        for i in range(5):
            self.assertAlmostEqual( g[i], true_values[i], delta=delta )
    
    
    def identical_strings( self, string, match=0, mismatch = -1, gap=-1, extension=-1 ):
        print("Testing '%s' against itself with parameters (%s, %s, %s, %s)" % (string, match, mismatch, gap, extension) )
        g = self.scores( string, string, match, mismatch, gap, extension )
        
        self.assertAlmostEqualGotoh( g, (match*len(string), len(string), 0,0,0), delta=1.0 )
    
    def replacement( self, string, char, match=0, mismatch = -1, gap=-1, extension=-1 ):
        print("Testing replacing '%s' in '%s'" % (char, string))
        
        n = string.count(char)
        replacement_string = string.replace(char, '_')
        print('%d of %d chars replaced' % (n, len(string)) )
        print("Replaced string:", replacement_string)
        g = self.scores( string, replacement_string, match, mismatch, gap, extension )
        print('scores returned:', g)
        self.assertAlmostEqualGotoh( g, (match*(len(string)-n)+mismatch*n, len(string)-n, n, 0,0), delta=0.5 )
    
    def test_scores(self):
        lorem = "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum"
        john_3_16_ar = "لانه هكذا احب الله العالم حتى بذل ابنه الوحيد لكي لا يهلك كل من يؤمن به بل تكون له الحياة الابدية"
        self.identical_strings(lorem)
        self.identical_strings(lorem, 2, -1, -5, -3)
        self.identical_strings(lorem, 4.2, -0.3, -5.2, 0.0)
        self.identical_strings(john_3_16_ar)
        self.identical_strings(john_3_16_ar, 2.3, -0.3, -0.2, -0.1)
        
        self.levenshtein( lorem, lorem.replace( 'o', 'r' ) )
        self.levenshtein( "AAA", "AABC" )
        
        self.replacement( lorem, 'r' )
        self.replacement( lorem, 'c', 4.2, -0.03, -0.1, -0.0003 )
            

if __name__ == "__main__":
    unittest.main()

