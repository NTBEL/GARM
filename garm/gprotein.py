# PySB components
from pysb import Model, Monomer, Parameter, Initial, Rule, Observable, Expression, Annotation, Compartment, ANY
# PySB macros
from pysb.macros import bind, bind_complex, catalyze, catalyze_complex, catalyze_state, degrade
# PySB alias components
from pysb.util import alias_model_components

def monomers():
    Monomer('G', ['b']) # G-protein (heterotrimer)
    Monomer('Ga',['bp','bgb', 'br']) # G_alpha subunit
    Monomer('Gb', ['bga','bgg']) # G_beta subunit
    Monomer('Gg', ['bgb']) # G_gamma subunit
    Monomer('GDP', ['b']) # GDP
    Monomer('GTP', ['b']) # GTP
    alias_model_components()
    return

def compartments(v_c=100):
    Parameter("V_C", v_extra)
    alias_model_components()
    Compartment('C', dimension=3, size=V_C) # C = Cytosol

def initials(g=1e4, gdp=200, gtp=500):
    Parameter('G_0', g) # G-protein trimer
    Paramter('GDP_0', gdp) # GDP
    Parameter('GTP_0', gtp) # GTP
    alias_model_components()
    Initial(G(b=None)**CM, G_0)
    Initial(GDP(b=None)**C, GDP_0)
    Initial(GTP(b=None)**C, GTP_0)
    return

def observables():
    alias_model_components()
    return

def ternary_complex(kd1=100, kf1=1e-3,
                    kd2=200, kf2=1e-4
                    kf3=1e-1, kr3=1,
                    kf4=1e-1, kr4=0.1,
                    kd5=100, kf5=0.01,
                    kd6=80, kf6=0.001):

    # L + R_I-G <--> L-R_I-G ; kd1, kf1
    # L + R_A-G <--> L-R_A-G ; kd2, kf2
    # L-R_I-G <--> L-R_A-G ; kf3, kr3
    # R_I-G <--> R_A-G ; kf4, kr4
    # R_I + G <---> R_I-G ; kd5, kf5
    # R_A + G <---> R_A-G; kd6, kf6
    Parameter('Kd_LRG_I', kd1)
    Parameter('kf_L_bind_R_IG', kf1)
    Parameter('Kd_LRG_A', kd2)
    Parameter('kf_L_bind_R_AG', kf2)
    Parameter('kf_LR_IG_to_A', kf3)
    Parameter('kr_LR_IG_to_A', kr3)
    Parameter('kf_R_IG_to_A', kf4)
    Parameter('kr_R_IG_to_A', kr4)
    Parameter('Kd_R_IG', kd5)
    Parameter('kf_G_bind_R_I', kf5)
    Parameter('Kd_R_AG', kd6)
    Parameter('kf_G_bind_R_A', kf6)
    alias_model_components()
    Expression('kr_L_bind_R_IG', Kd_LRG_I*kf_L_bind_R_IG)
    Expression('kr_L_bind_R_AG', Kd_LRG_A*kf_L_bind_R_AG)
    Expression('kr_G_bind_R_I', Kd_R_IG*kf_G_bind_R_I)
    Expression('kr_G_bind_R_A', Kd_R_AG*kf_G_bind_R_A)
    alias_model_components()

    L_free = L(b=None)**EXTRA
    R_I = R(bl=1, bg=None, state='I')**CM # Alias the inactive receptor
    R_A = R(bl=1, bg=None, state='A')**CM # Alias the active receptor
    LR_I = L(b=1)**EXTRA % R(bl=1, bg=None, state='I')**CM # Alias the inactive receptor complex
    LR_A = L(b=1)**EXTRA % R(bl=1, bg=None, state='A')**CM # Alias the active receptor complex
    R_I_G = R(bl=None, bg=2, state='I')**CM % G(b=2)**CM # Alias the inactive receptor-Gprotein complex
    R_A_G = R(bl=None, bg=2, state='A')**CM % G(b=2)**CM # Alias the inactive receptor-Gprotein complex
    LR_I_G = L(b=1)**EXTRA % R(bl=1, bg=2, state='I')**CM % G(b=2)**CM # Alias the inactive ligand-receptor-Gprotein complex
    LR_A_G = L(b=1)**EXTRA % R(bl=1, bg=2, state='A')**CM % G(b=2)**CM # Alias the active ligand-receptor-Gprotein complex
    # L + R_I-G <--> L-R_I-G ; kd1, kf1
    Rule(L_free + R_I_G | LR_I_G, kf_L_bind_R_IG, kr_L_bind_R_IG)
    # L + R_A-G <--> L-R_A-G ; kd2, kf2
    Rule(L_free + R_A_G | LR_A_G, kf_L_bind_R_AG, kr_L_bind_R_AG)
    # L-R_I-G <--> L-R_A-G ; kf3, kr3
    Rule(LR_I_G | LR_A_G, kf_LR_IG_to_A, kr_LR_IG_to_A)
    # R_I-G <--> R_A-G ; kf4, kr4
    Rule(R_I_G | R_A_G, kf_R_IG_to_A, kr_R_IG_to_A)
    # R_I + G <---> R_I-G ; kd5, kf5
    Rule(R_I + G(b=None)**CM, kf_G_bind_R_I, kr_G_bind_R_I)
    # R_A + G <---> R_A-G; kd6, kf6
    Rule(R_A + G(b=None)**CM, kf_G_bind_R_A, kr_G_bind_R_A)
    return
