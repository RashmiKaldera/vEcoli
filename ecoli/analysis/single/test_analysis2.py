import os
import numpy as np
from typing import Any

from duckdb import DuckDBPyConnection
import altair as alt
import pandas as pd


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
    SELECT time, listeners__growth_limits__ppgpp_conc
    FROM ({history_sql})
    ORDER BY time ASC
    """

    output_df = conn.sql(query).df()
    # Convert time from seconds to minutes
    output_df["Time (min)"] = output_df["time"] / 60

    # Create Altair line chart
    chart = (
        alt.Chart(output_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__growth_limits__ppgpp_conc:Q",
                title="ppGpp Concentration (µM)",
            ),
        )
        .properties(
            title="ppGpp Concentration Over Time",
            width=600,
            height=400,
        )
    )
    output_df.to_csv(os.path.join(outdir, "ppgpp_concentration.csv"), index=False)
    html_path = os.path.join(outdir, "ppgpp_concentration.html")
    chart.save(html_path)

    query2 = f"""
        SELECT time, bulk
        FROM ({history_sql})
        ORDER BY time ASC
        """
    output_df2 = conn.sql(query2).df()
    time_minutes = output_df2["time"].to_numpy() / 60
    bulk_matrix = np.stack(output_df2["bulk"].values).astype(int)
    bulk_df = pd.DataFrame(bulk_matrix)
    bulk_df.insert(0, "Time (min)", time_minutes)
    bulk_df.to_csv(os.path.join(outdir, "bulk_matrix.csv"), index=False)

    # deoxyribose
    deoxyribose_col = bulk_df.columns[70]
    # Create a DataFrame with just Time and Deoxyribose
    deoxyribose_df = bulk_df[["Time (min)", deoxyribose_col]].copy()
    deoxyribose_df.rename(columns={deoxyribose_col: "Deoxyribose"}, inplace=True)

    # # Step 5: Create new DataFrame for plotting

    chart2 = (
        alt.Chart(deoxyribose_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Deoxyribose:Q", title="Deoxyribose (counts)"),
        )
        .properties(title="Deoxyribose over Time")
    )
    chart2.save(os.path.join(outdir, "deoxyribose_concentration.html"))

    # Select time and the 3 ammonium-related columns
    ammonium_cols = [655, 656, 657]
    selected_cols = ["Time (min)"] + [bulk_df.columns[i] for i in ammonium_cols]
    ammonium_df = bulk_df[selected_cols].copy()

    # Rename for clarity
    ammonium_df.columns = ["Time (min)", "Ammonium (C)", "Ammonium (E)", "Ammonium (P)"]

    # Melt for plotting
    melted_ammonium_df = pd.melt(
        ammonium_df, id_vars=["Time (min)"], var_name="Compartment", value_name="Counts"
    )

    # Plot
    chart3 = (
        alt.Chart(melted_ammonium_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Counts:Q", title="Ammonium (counts)"),
            color=alt.Color("Compartment:N", scale=alt.Scale(range=COLORS)),
        )
        .properties(title="Ammonium Counts Over Time")
    )

    # Save or display
    chart3.save(os.path.join(outdir, "ammonium_counts.html"))

    # Select time and the 3 ATP-related columns
    atp_cols = [797, 798, 799]
    selected_atp_cols = ["Time (min)"] + [bulk_df.columns[i] for i in atp_cols]
    atp_df = bulk_df[selected_atp_cols].copy()

    # Rename for clarity
    atp_df.columns = ["Time (min)", "ATP (C)", "ATP (P)", "ATP (E)"]

    # Melt for plotting
    melted_atp_df = pd.melt(
        atp_df, id_vars=["Time (min)"], var_name="Compartment", value_name="Counts"
    )

    # Plot
    chart4 = (
        alt.Chart(melted_atp_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Counts:Q", title="ATP (counts)"),
            color=alt.Color("Compartment:N", scale=alt.Scale(range=COLORS)),
        )
        .properties(title="ATP Counts Over Time")
    )

    # Save or display
    chart4.save(os.path.join(outdir, "ATP_counts.html"))

    # Select time and the 3 Glucose-related columns
    glucose_cols = [6105, 6106]
    selected_glucose_cols = ["Time (min)"] + [bulk_df.columns[i] for i in glucose_cols]
    glucose_df = bulk_df[selected_glucose_cols].copy()

    # Rename for clarity
    glucose_df.columns = ["Time (min)", "D-Glucose (C)", "D-Glucose (P)"]

    # Melt for plotting
    melted_glucose_df = pd.melt(
        glucose_df, id_vars=["Time (min)"], var_name="Compartment", value_name="Counts"
    )

    # Plot
    chart5 = (
        alt.Chart(melted_glucose_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Counts:Q", title="Glucose (counts)"),
            color=alt.Color("Compartment:N", scale=alt.Scale(range=COLORS)),
        )
        .properties(title="Glucose Counts Over Time")
    )

    # Save or display
    chart5.save(os.path.join(outdir, "Glucose_counts.html"))

    # glycogen
    glycogen_col = bulk_df.columns[9945]
    # Create a DataFrame with just Time and glycogen
    glycogen_df = bulk_df[["Time (min)", glycogen_col]].copy()
    glycogen_df.rename(columns={glycogen_col: "Glycogen"}, inplace=True)

    # # Step 5: Create new DataFrame for plotting

    chart6 = (
        alt.Chart(glycogen_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Glycogen:Q", title="Glycogen (counts)"),
        )
        .properties(title="Glycogen over Time")
    )
    chart6.save(os.path.join(outdir, "glycogen_counts.html"))

    query3 = f"""
        SELECT time, listeners__rna_counts__mRNA_counts,listeners__rna_counts__full_mRNA_counts,listeners__rna_counts__partial_mRNA_counts,listeners__rna_counts__mRNA_cistron_counts,listeners__rna_counts__full_mRNA_cistron_counts,listeners__rna_counts__partial_mRNA_cistron_counts,listeners__rna_counts__partial_rRNA_counts,listeners__rna_counts__partial_rRNA_cistron_counts
        FROM ({history_sql})
        ORDER BY time ASC
        """

    output_df3 = conn.sql(query3).df()
    # Convert time from seconds to minutes
    output_df3["Time (min)"] = output_df3["time"] / 60

    mrna_matrix = np.stack(
        output_df3["listeners__rna_counts__mRNA_counts"].values
    ).astype(int)
    full_matrix = np.stack(
        output_df3["listeners__rna_counts__full_mRNA_counts"].values
    ).astype(int)
    partial_matrix = np.stack(
        output_df3["listeners__rna_counts__partial_mRNA_counts"].values
    ).astype(int)

    rna_df = pd.DataFrame(
        {
            "Time (min)": output_df3["Time (min)"],
            "mRNA": mrna_matrix[:, 0],
            "full mRNA": full_matrix[:, 0],
            "partial mRNA": partial_matrix[:, 0],
        }
    )

    melted_rna_df = pd.melt(
        rna_df, id_vars=["Time (min)"], var_name="RNA Type", value_name="Counts"
    )

    # Plot
    chart_rna = (
        alt.Chart(melted_rna_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Counts:Q", title="RNA Counts"),
            color=alt.Color(
                "RNA Type:N", title="RNA Type", scale=alt.Scale(range=COLORS)
            ),
            # Optional color customization
        )
        .properties(title="RNA Counts Over Time", width=700, height=400)
    )

    chart_rna.save(os.path.join(outdir, "rna_counts.html"))

    query4 = f"""
        SELECT time, listeners__unique_molecule_counts__active_RNAP
        FROM ({history_sql})
        ORDER BY time ASC
        """

    output_df4 = conn.sql(query4).df()
    # Convert time from seconds to minutes
    output_df4["Time (min)"] = output_df4["time"] / 60

    # Create Altair line chart
    chart_rnap = (
        alt.Chart(output_df4)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__active_RNAP:Q",
                title="Active RNAP count",
            ),
        )
        .properties(
            title="Active RNAP count Over Time",
            width=600,
            height=400,
        )
    )
    output_df4.to_csv(os.path.join(outdir, "active_RNAP.csv"), index=False)
    chart_rnap.save(os.path.join(outdir, "active_RNAP.html"))

    query5 = f"""
            SELECT time, listeners__unique_molecule_counts__RNA
            FROM ({history_sql})
            ORDER BY time ASC
            """

    output_df5 = conn.sql(query5).df()
    # Convert time from seconds to minutes
    output_df5["Time (min)"] = output_df5["time"] / 60

    # Create Altair line chart
    chart_rnau = (
        alt.Chart(output_df5)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("listeners__unique_molecule_counts__RNA:Q", title="RNA count"),
        )
        .properties(
            title="RNA count Over Time",
            width=600,
            height=400,
        )
    )
    output_df5.to_csv(os.path.join(outdir, "rna_unique.csv"), index=False)
    chart_rnau.save(os.path.join(outdir, "RNAunique.html"))

    query6 = f"""
                SELECT time, listeners__growth_limits__aa_conc
                FROM ({history_sql})
                ORDER BY time ASC
                """

    output_df6 = conn.sql(query6).df()
    # Convert time from seconds to minutes
    output_df6["Time (min)"] = output_df6["time"] / 60

    # Create Altair line chart
    chart_aa = (
        alt.Chart(output_df6)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__growth_limits__aa_conc:Q", title="Amino acid concentration"
            ),
        )
        .properties(
            title="Amino Acid concentration Over Time",
            width=600,
            height=400,
        )
    )
    output_df6.to_csv(os.path.join(outdir, "Aminoacidconc.csv"), index=False)
    chart_aa.save(os.path.join(outdir, "Aminoacidconc.html"))

    query7 = f"""
                    SELECT time, listeners__unique_molecule_counts__oriC
                    FROM ({history_sql})
                    ORDER BY time ASC
                    """

    output_df7 = conn.sql(query7).df()
    # Convert time from seconds to minutes
    output_df7["Time (min)"] = output_df7["time"] / 60

    # Create Altair line chart
    chart_oric = (
        alt.Chart(output_df7)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("listeners__unique_molecule_counts__oriC:Q", title="Count of oriC"),
        )
        .properties(
            title="oriC count Over Time",
            width=600,
            height=400,
        )
    )

    chart_oric.save(os.path.join(outdir, "oric.html"))
