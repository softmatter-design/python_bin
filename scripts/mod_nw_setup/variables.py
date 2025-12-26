###########
# basic_cond = [ver_cognac, blank_udf, base_udf, core]
ver_cognac = ''
blank_udf = ''
base_udf = ''
core = 1

# nw_cond = [nw_model, strand, n_strand, n_segments, n_cell, n_sc, l_bond, c_n]
nw_model = ''
calc = ''
strand = ''
n_strand = 1
n_segments = 1
strand_type = ''
n_cell = 1
n_sc = 1
l_bond = 1.0
c_n = 1.0

# rnd_cond = [restart, cond_top, histgram_bins]
restart = ''
cond_top = []
histgram_bins = 0

# sim_cond = [entanglement, multi_org, density_org, shrinkage, expand, step_press, press_time, step_rfc, step_rfc_time, eqn_repeat, eqn_time, greenkubo, greenkubo_repeat, greenkubo_time, calc]
entanglement = ''
multi_org = 0
density_org = 0
shrinkage = 0
expand = 1.0
step_press = []
press_time = []
step_rfc = []
step_rfc_time = []
eqn_repeat = 1
eqn_time = []
greenkubo = ''
greenkubo_repeat = 1
greenkubo_time = []

exchange1=''
exchange1_repeat = 1
exchange1_target=''
exchange1_cond=[]
exchange1_post_time = []
exchange1_eqn_repeat = 1
exchange1_eqn_time = []

exchange2=''
exchange2_repeat = 1
exchange2_target=''
exchange2_cond=[]
exchange2_post_time = []
exchange2_eqn_repeat = 1
exchange2_eqn_time = []


# mod_cond = [n_chains, n_beads_unit, e2e, org_unitcell]
n_chains = 1
n_jp = 1
n_beads_unit = 1
e2e = 1.0
org_unitcell = 1.0

# sim_cond2 = [multi_mod, density_mod]
multi_mod = 1
density_mod = 1.0

# target_cond = [system, unit_cell, total_net_atom, nu]
system = 1.0
unit_cell = 1.0
total_net_atom = 1
nu = 1.0
d_jp = 1

target_dir = ''
target_name = ''


# Cognac用の名称設定
nw_name = "Network"
atom_name = ["JP_A", "Strand_A", "Dis_D", "Dis_S"]
bond_name = ["bond_JP-Chn", "bond_Strand", "bond_DD", "bond_DS"]
angle_name = ["angle_AAA"]
site_name = ["site_JP", "site_Strand", "site_Dis_D", "site_Dis_S"]
pair_name = ["site_JP-site_JP", "site_JP-site_Strand", "site_JP-site_Dis_D", "site_JP-site_Dis_S",
            "site_Strand-site_Strand", "site_Strand-site_Dis_D", "site_Strand-site_Dis_S",
            "site_Dis_D-site_Dis_D", "site_Dis_D-site_Dis_S",
            "site_Dis_S-site_Dis_S"
            ]
site_pair_name = [ 
                ["site_JP", "site_JP"], 
                ["site_JP","site_Strand"], 
                ["site_JP", "site_Dis_D"], 
                ["site_JP", "site_Dis_S"], 
                ["site_Strand", "site_Strand"],
                ["site_Strand", "site_Dis_D"], 
                ["site_Strand", "site_Dis_S"], 
                ["site_Dis_D", "site_Dis_D"], 
                ["site_Dis_D", "site_Dis_S"], 
                ["site_Dis_S", "site_Dis_S"]
                ]
# Set [R0, K]
harmonic = [0.97, 1000]
# [Potential_Type, theta0, K]		
angle = ['Theta2', 74, 10.0]
# [Cutoff, Scale_1_4_Pair, sigma, epsilon, range]
lj_cond = [2**(1/6), 1.0, 1.0, 1.0, 1.0]	

# Condition for Dissociative Bond
n_spacer = 2