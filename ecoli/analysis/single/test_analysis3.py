from ecoli.library.schema import bulk_name_to_idx
import pandas as pd

# %%
# Loading Simdata to obtain the bulk id labels
from ecoli.library.sim_data import LoadSimData

sim_data_default = "out/plasmid/parca/kb/simData.cPickle"
sim_data = LoadSimData(sim_data_default).sim_data

bulk_molecule_ids = sim_data.internal_state.bulk_molecules.bulk_data[
    "id"
].tolist()  # Model common name with compartments
# %%


df1 = pd.DataFrame(bulk_molecule_ids, columns=["bulk_id"])
df1.to_csv("bulk_molecule_ids1.csv", index=False)
# %%
dntp_ids = bulk_name_to_idx(
    ["DATP[c]", "DCTP[c]", "DGTP[c]", "TTP[c]"], bulk_molecule_ids
)
polymerized_dntp_ids = bulk_name_to_idx(
    [
        "polymerized_DATP[c]",
        "polymerized_DCTP[c]",
        "polymerized_DGTP[c]",
        "polymerized_TTP[c]",
    ],
    bulk_molecule_ids,
)
ntp_ids = bulk_name_to_idx(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"], bulk_molecule_ids)
polymerized_ntp_ids = bulk_name_to_idx(
    [
        "polymerized_ATP[c]",
        "polymerized_CTP[c]",
        "polymerized_GTP[c]",
        "polymerized_UTP[c]",
    ],
    bulk_molecule_ids,
)
