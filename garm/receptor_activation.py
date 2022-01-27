# PySB components
from pysb import Model, Monomer, Parameter, Initial, Rule, Observable, Expression, Annotation, Compartment, ANY
# PySB macros
from pysb.macros import bind, bind_complex, catalyze, catalyze_complex, catalyze_state, degrade
# PySB alias components
from pysb.util import alias_model_components


def monomers():
    # bl = binding site for Ligand
    # bg = binding site for G-protein
    Monomer('R',['bl', 'bg','state'], {'state': ['I','A']}) # Receptor
    Monomer('L', ['b']) # Ligand
    alias_model_components()
    return

def compartments(v_extra=1000, v_cm=100):
    Parameter("V_EXTRA", v_extra)
    Parameter("V_CM", v_cm)
    alias_model_components()
    Compartment('EXTRA', dimension=3, size=V_EXTRA)
    # Cell Membrane
    Compartment('CM', dimension=2, parent=EXTRA, size=V_CM)

def initials(l=1e4, r=200):
    Parameter('L_0', l) # Ligand
    Paramter('R_0', r) # Receptor
    alias_model_components()
    Initial(L(b=None)**EXTRA, L_0)
    Initial(R(bl=None, bg=None, state='I')**CM, R_0)
    return

def observables():
    alias_model_components()
    Observable('L_free',L(b=None)**EXTRA) # Free ligand
    Observable('LR', L(b=1)**EXTRA % R(bl=1, bg=None)**CM) # total ligand-receptor complex
    Observable('LR_I', L(b=1)**EXTRA % R(bl=1, bg=None, state='I')**CM) # ligand-inactive receptor complex
    Observable('LR_A', L(b=1)**EXTRA % R(bl=1, bg=None, state='A')**CM) # ligand-active receptor complex
    Observable('R_I', R(state='I')**CM) # total inactive receptor
    Observable('R_A', R(state='A')**CM) # total active receptor
    return


def single_state_activation(kd=100, kf=1e-3):
    # L + R_I <--> LR_A
    Parameter('Kd_LR', kd)
    Parameter('kf_L_bind_R', kf)
    alias_model_components()
    Expression('kr_L_bind_R', Kd_LR*kf_L_bind_R)
    alias_model_components()
    LR_A = L(b=1)**EXTRA % R(bl=1, bg=None, state='A')**CM # Alias the active receptor complex
    Rule(L(b=None)**EXTRA + R(bl=None, bg=None, state='I') | LR_A, kf_L_bind_R, kr_L_bind_R)
    return

def two_state_activation(kd = 100, kf=1e-3, kfa=1e-1, kra=1):
    # L + R_I <--> LR_I <--> LR_A
    Parameter('Kd_LR', kd)
    Parameter('kf_L_bind_R', kf)
    Parameter('kf_R_I_to_A', kfa)
    Parameter('kr_R_I_to_A', kra)
    alias_model_components()
    Expression('kr_L_bind_R', Kd_LR*kf_L_bind_R)
    alias_model_components()
    LR_I = L(b=1)**EXTRA % R(bl=1, bg=None, state='I')**CM # Alias the inactive receptor complex
    LR_A = L(b=1)**EXTRA % R(bl=1, bg=None, state='A')**CM # Alias the active receptor complex
    Rule(L(b=None)**EXTRA + R(bl=None, bg=None, state='I') | LR_I, kf_L_bind_R, kr_L_bind_R)
    Rule(LR_I | LR_A, kf_R_I_to_A, kr_R_I_to_A)
    return

def reversible_two_state_activation(kd_i=100, kf_i=1e-3,
                                    kd_a=200, kf_a=1e-4
                                    kfa=1e-1, kra=1,
                                    kfia=1e-1, kria=0.1):

    # L + R_I <--> LR_I
    # L + R_A <--> LR_A
    # LR_I <--> LR_A
    # R_I <--> R_A
    Parameter('Kd_LR_I', kd_i)
    Parameter('kf_L_bind_R_I', kf_i)
    Parameter('Kd_LR_A', kd_a)
    Parameter('kf_L_bind_R_A', kf_a)
    Parameter('kf_LR_I_to_A', kfa)
    Parameter('kr_LR_I_to_A', kra)
    Parameter('kf_R_I_to_A', kfia)
    Parameter('kr_R_I_to_A', kria)
    alias_model_components()
    Expression('kr_L_bind_R_I', Kd_LR_I*kf_L_bind_R_I)
    Expression('kr_L_bind_R_A', Kd_LR_A*kf_L_bind_R_A)
    alias_model_components()
    R_I = R(bl=1, bg=None, state='I')**CM # Alias the inactive receptor
    R_A = R(bl=1, bg=None, state='A')**CM # Alias the active receptor
    LR_I = L(b=1)**EXTRA % R(bl=1, bg=None, state='I')**CM # Alias the inactive receptor complex
    LR_A = L(b=1)**EXTRA % R(bl=1, bg=None, state='A')**CM # Alias the active receptor complex
    Rule(L(b=None)**EXTRA + R(bl=None, bg=None, state='I') | LR_I, kf_L_bind_R_I, kr_L_bind_R_I)
    Rule(L(b=None)**EXTRA + R(bl=None, bg=None, state='A') | LR_A, kf_L_bind_R_A, kr_L_bind_R_A)
    Rule(LR_I | LR_A, kf_R_I_to_A, kr_R_I_to_A)
    Rule(R_I | R_A, kf_R_I_to_A, kr_R_I_to_A)
    return
