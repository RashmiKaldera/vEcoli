import os
from typing import Any

from duckdb import DuckDBPyConnection


COLORS_256 = [  # From colorbrewer2.org, qualitative 8-class set 1
    [228, 26, 28],
    [55, 126, 184],
    [77, 175, 74],
    [152, 78, 163],
    [255, 127, 0],
    [255, 255, 51],
    [166, 86, 40],
    [247, 129, 191],
]

COLORS = ["#%02x%02x%02x" % (color[0], color[1], color[2]) for color in COLORS_256]


def plot(
    params: dict[str, Any],
    conn: DuckDBPyConnection,
    history_sql: str,
    config_sql: str,
    success_sql: str,
    sim_data_paths: dict[str, dict[int, str]],
    validation_data_paths: list[str],
    outdir: str,
    variant_metadata: dict[str, dict[int, Any]],
    variant_names: dict[str, str],
):
    query = f"""
    SELECT * FROM ({history_sql})
    ORDER BY time ASC;
    """

    output_df = conn.sql(query).df()
    # bulk_matrix = np.stack(output_df["bulk"].values).astype(int)
    # np.savetxt(os.path.join(outdir, "bulk_matrix.txt"), bulk_matrix)
    with open(os.path.join(outdir, "output_df_columns.txt"), "w") as f:
        f.write("\n".join(output_df.columns))

    # with open(os.path.join(outdir, "history_sql.txt"), "w") as f:
    # f.write(history_sql)
    #
    #
    # query_dict = {
    #     "experiment_id": "plasmidtest",
    #     "variant": 0,
    #     "lineage_seed": 0,
    #     "generation": 1
    # }
    # query = f'''
    # SELECT *,time FROM read_parquet("out/{query_dict["experiment_id"]}/history/*/*/*/*/*/*.pq", hive_partitioning=true)
    # WHERE variant={query_dict['variant']}
    # AND lineage_seed={query_dict['lineage_seed']}
    # AND generation={query_dict['generation']}
    # ORDER BY time
    # '''
