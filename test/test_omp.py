from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import unittest as ut
import numpy.testing as nt 

from wrf import (omp_set_num_threads, omp_get_num_threads, 
                 omp_get_max_threads, omp_get_thread_num,
                 omp_get_num_procs, omp_in_parallel, 
                 omp_set_dynamic, omp_get_dynamic, omp_set_nested,
                 omp_get_nested, omp_set_schedule, 
                 omp_get_schedule, omp_get_thread_limit, 
                 omp_set_max_active_levels, 
                 omp_get_max_active_levels, omp_get_level,
                 omp_get_ancestor_thread_num, omp_get_team_size,
                 omp_get_active_level, omp_in_final,
                 omp_init_lock, omp_init_nest_lock,
                 omp_destroy_lock, omp_destroy_nest_lock,
                 omp_set_lock, omp_set_nest_lock,
                 omp_unset_lock, omp_unset_nest_lock,
                 omp_test_lock, omp_test_nest_lock,
                 omp_get_wtime, omp_get_wtick)
from wrf import Constants
            

class OmpTest(ut.TestCase):
    longMessage = True
    
    def test_locks(self):
        l = omp_init_lock()
        omp_set_lock(l)
        omp_unset_lock(l)
        omp_test_lock(l)
        omp_destroy_lock(l)
        
        nl = omp_init_nest_lock()
        omp_set_nest_lock(nl)
        omp_unset_nest_lock(nl)
        omp_test_nest_lock(nl)
        omp_destroy_nest_lock(nl)
    
    
    def test_thread_set(self):
        omp_set_num_threads(4)
        max_threads = omp_get_max_threads()
        self.assertEqual(max_threads, 4)
        
        num_threads = omp_get_num_threads()
        self.assertEqual(num_threads, 1) # Always 1 outside of parallel region
        
        thread_num = omp_get_thread_num()
        self.assertEqual(thread_num, 0) # Always 0 outside of parallel region
        num_procs = omp_get_num_procs()
        in_parallel = omp_in_parallel()
        self.assertFalse(in_parallel) # Always False outside of parallel region
        
        limit = omp_get_thread_limit()
        
    
    def test_dynamic(self):
        omp_set_dynamic(True)
        dynamic = omp_get_dynamic()
        self.assertTrue(dynamic)
        
        omp_set_dynamic(False)
        dynamic = omp_get_dynamic()
        self.assertFalse(dynamic)
    
    def test_nested(self):
        omp_set_nested(True)
        nested = omp_get_nested()
        self.assertTrue(nested)
        
        omp_set_nested(False)
        nested = omp_get_nested()
        self.assertFalse(nested)
        
             
    def test_schedule(self):
        omp_set_schedule(Constants.OMP_SCHED_STATIC, 100000)
        kind, modifier = omp_get_schedule()
        self.assertEqual(kind, Constants.OMP_SCHED_STATIC)
        self.assertEqual(modifier, 100000)
        
        omp_set_schedule(Constants.OMP_SCHED_DYNAMIC, 10000)
        kind, modifier = omp_get_schedule()
        self.assertEqual(kind, Constants.OMP_SCHED_DYNAMIC)
        self.assertEqual(modifier, 10000)
        
        omp_set_schedule(Constants.OMP_SCHED_GUIDED, 100) 
        kind, modifier = omp_get_schedule()
        self.assertEqual(kind, Constants.OMP_SCHED_GUIDED)
        self.assertEqual(modifier, 100)
        
        omp_set_schedule(Constants.OMP_SCHED_AUTO, 10) 
        kind, modifier = omp_get_schedule()
        self.assertEqual(kind, Constants.OMP_SCHED_AUTO)
        self.assertNotEqual(modifier, 10) # The modifier argument is ignored,
                                          # so it will be set to the previous
                                          # value of 100.
        
    
    def test_team_level(self):
        omp_set_max_active_levels(10)
        active_levels = omp_get_max_active_levels()
        self.assertEqual(active_levels, 10)
        
        level = omp_get_level()
        ancestor_thread = omp_get_ancestor_thread_num(level)
        team_size = omp_get_team_size(level)
        active_level = omp_get_active_level()
        in_final = omp_in_final()
        
    
    def test_time(self):
        wtime = omp_get_wtime()
        wtick = omp_get_wtick()

if __name__ == "__main__":
    ut.main()
    