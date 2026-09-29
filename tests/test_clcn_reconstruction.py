"""Focused discriminator test: short insertion versus deletion-like left join."""
import importlib.util
import random
from pathlib import Path
from types import SimpleNamespace
import unittest

spec=importlib.util.spec_from_file_location('reconstruct',Path(__file__).parents[1]/'scripts/clcnkb/clcn_reconstruct.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class ReconstructionTest(unittest.TestCase):
    def test_both_junctions_deduplicate_and_filter_quality(self):
        rng=random.Random(91)
        seq=lambda n:''.join(rng.choice('ACGT') for _ in range(n))
        left,right,copied=seq(100),seq(100),seq(79)
        ins=left[-35:]+copied+right[:35]
        deletion_like=left[-35:]+copied+seq(35)
        discriminator=left[-10:]+copied+right[:10]
        self.assertIn(discriminator,ins)
        self.assertNotIn(discriminator,deletion_like)
        probe=m.canon(left[-15:]+copied[:16])
        a=SimpleNamespace(query_sequence=ins,query_qualities=[30]*len(ins))
        b=SimpleNamespace(query_sequence=m.rc(ins),query_qualities=[30]*len(ins))
        # Same template's two overlapping mates must contribute only once.
        hits=set(m.eligible_kmers(a,31))|set(m.eligible_kmers(b,31))
        self.assertIn(probe,hits)
        self.assertEqual(sum(x==probe for x in hits),1)
        a.query_qualities[30]=10
        self.assertNotIn(probe,set(m.eligible_kmers(a,31)))
        self.assertEqual(m.canon(ins[:51]),m.canon(m.rc(ins[:51])))

if __name__=='__main__':unittest.main()
