import numpy as np

class Contig:
    def __init__(self, con, length):
        self.name = con
        self.length = 0
        self.coverage = np.array([])
        self.abundance = np.array([])
    
    def get_mean(self):
        return np.mean(self.coverage)
    def get_std(self):
        return np.std(self.coverage)

class Window:
    def __init__(self, con, start, end, depth):
        self.contig = [con]
        self.start = start
        self.end = end
        self.depth = depth
        self.length = end-start+1
        self.coverage = np.array([])
        self.abundance = np.array([])
        self.cov_mean = 0
        self.abund_mean = 0
        self.aligns = []
        self.profile = np.array([])
        self.align_locs = []
    def get_abund_mean(self):
        return np.mean(self.abundance)
    def get_abund_std(self):
        return np.std(self.abundance)
    def get_cov_mean(self):
        return np.mean(self.coverage)
    def get_cov_std(self):
        return np.std(self.coverage)

class subContig:
    def __init__(self, con, start, end, ref):
        self.contig = con
        self.start = start
        self.end = end
        self.ref = ref
        self.length = end-start+1
        self.coverage = np.array([])
        self.abundance = np.array([])
        self.cov_mean = 0
        self.abund_mean = 0
        self.label = 0
        self.probs = np.array([])
        self.priors = np.array([])
    def get_abund_mean(self):
        return np.mean(self.abundance)
    def get_abund_std(self):
        return np.std(self.abundance)
    def get_cov_mean(self):
        return np.mean(self.coverage)
    def get_cov_std(self):
        return np.std(self.coverage)
