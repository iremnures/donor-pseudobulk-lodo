import pandas as pd

cells = pd.DataFrame(
    [
        ["D1", "D1_control", "control", "Monocyte", 2, 4],
        ["D1", "D1_control", "control", "Monocyte", 0, 3],
        ["D1", "D1_control", "control", "B cell",   1, 8],
        ["D1", "D1_treated", "treated", "Monocyte", 6, 4],
        ["D1", "D1_treated", "treated", "Monocyte", 4, 2],
        ["D1", "D1_treated", "treated", "B cell",   2, 7],
        ["D2", "D2_control", "control", "Monocyte", 3, 5],
        ["D2", "D2_control", "control", "B cell",   0, 9],
        ["D2", "D2_treated", "treated", "Monocyte", 8, 3],
        ["D2", "D2_treated", "treated", "B cell",   1, 8],
    ],
    columns=[
        "donor_id", "sample_id", "condition",
        "cell_type", "GENE_A", "GENE_B"
    ]
)

keys = ["donor_id", "sample_id", "condition", "cell_type"]

pseudobulk = cells.groupby(keys, sort=True)[
    ["GENE_A", "GENE_B"]
].sum()

cell_numbers = cells.groupby(keys, sort=True).size().rename("n_cells")
summary = pseudobulk.join(cell_numbers)

assert pseudobulk["GENE_A"].sum() == cells["GENE_A"].sum()
assert pseudobulk["GENE_B"].sum() == cells["GENE_B"].sum()

print("--- Cell-level raw counts ---")
print(cells.to_string(index=False))

print("\n--- Pseudobulk raw counts and cell coverage ---")
print(summary.to_string())

print("\nIndependent donors:", cells["donor_id"].nunique())
print("Biological samples:", cells["sample_id"].nunique())
print("Measured cells:", len(cells))

summary.to_csv("results/tables/toy_pseudobulk_summary.csv")
print("\nToy summary saved. No statistical testing performed.")
