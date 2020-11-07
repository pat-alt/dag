## DAG module

class Dag():
    def __init__(self, scm, outcome, exposure):
        self.scm = scm # list of dictionaries
        self.outcome = outcome # outcome variable
        self.exposure = exposure # exposure variable

    def find_paths(self):
        pass

    def adjustment_sets(self):
        pass

    def is_valid_adjustment_set(self, proposed_set):
        pass
