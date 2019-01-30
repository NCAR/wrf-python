from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from threading import Thread
try:
    from Queue import Queue
except ImportError:
    from queue import Queue
from collections import OrderedDict

import unittest as ut
import numpy.testing as nt

from wrf.cache import cache_item, get_cached_item, _get_cache
from wrf.config import get_cache_size


class TestThread(Thread):
    def __init__(self, num, q):
        self.num = num
        self.q = q
        super(TestThread, self).__init__()

    def run(self):
        for i in range(get_cache_size() + 10):
            key = "A" + str(i)
            cache_item(key, "test", i * self.num)

            item = get_cached_item(key, "test")

            if item != i * self.num:
                raise RuntimeError("cache is bogus")

        cache = OrderedDict(_get_cache())

        self.q.put(cache)


class CacheTest(ut.TestCase):
    longMessage = True

    def test_thread_local(self):
        q1 = Queue()
        q2 = Queue()
        thread1 = TestThread(2, q1)
        thread2 = TestThread(40, q2)

        thread1.start()
        thread2.start()

        result1 = q1.get(True, 1)
        result2 = q2.get(True, 1)

        thread1.join()
        thread2.join()

        print(result1)
        print(result2)

        # Result 1 and 2 shoudl be different
        self.assertNotEqual(result1, result2)

        # This thread should have no cache
        self.assertIsNone(_get_cache())


if __name__ == "__main__":
    ut.main()
