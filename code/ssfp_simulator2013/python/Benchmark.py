import time

class Timer:    
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

def test():
    for n in range(100000):
        a = 10 * 10

if __name__ == '__main__':
    import MRIpy
    import MRI
    import math

    try:
        with Timer() as t:
            MRIpy.SpectrumTest(True)
            #MRIpy.ImageTest(True)
    finally:
        print t.interval