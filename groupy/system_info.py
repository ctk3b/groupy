class System_info():
    def __init__(self, t_components=None):
        self.n_components = 0
        self.n_molecules = 0
        self.components = list()
        if t_components:
            self.add_components(t_components)

    def add_components(self, t_components):
        for t_component in t_components:
            assert len(t_component) == 3
            self.components.append(t_component)
            self.n_components += 1
            self.n_molecules += t_component[0]

    def print_system_info(self, filename="system-composition.txt"):
        f = open(filename, 'w')
        f.write('SYSTEM COMPOSITION\n\n')
        f.write('Species     Molecules     Atoms per molecule\n')
        f.write('-------     ---------     ------------------\n')
        for i, component in enumerate(self.components):
            f.write('%-12s%-14d%-d\n' % (component[2], component[0], component[1]))
        f.write('\n')
        f.write('%d total molecules' % self.n_molecules)
