from modules.auxiliary_modules.cDNA_generation import cDNA_generator

class primer_design:
    def __init__(self, u, l, o="./Outputs", b=59, v=2, b1 = 58.5, v1=2,b2=71, v2=2, dg_min=-4.5, dg_max=4.5, dg_min_FR = -5.5, dg_max_FR = 5.5,dg_min_PP = -7.5, dg_max_PP = 7.5,temp=37,**kwargs):
        self.cDNA = cDNA_generator(u, l).fetch_sequence()
        self.output_dir = o
        self.preserved_sRNA = u
        self.preserved_stemloop = l
        self.temp = int(temp)
        self.tm_min = b - v
        self.tm_max = b + v
        self.tm_min_probe_1 = b1 - v1
        self.tm_max_probe_1 = b1 + v1
        self.tm_min_probe_2 = b2 - v2
        self.tm_max_probe_2 = b2 + v2
        self.dg_min = dg_min
        self.dg_max = dg_max
        self.dg_min_FR = dg_min_FR,
        self.dg_max_FR = dg_max_FR,
        self.dg_min_PP = dg_min_PP,
        self.dg_max_PP = dg_max_PP