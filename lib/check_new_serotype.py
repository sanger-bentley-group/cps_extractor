import pandas as pd


class NewSerotype:
    def __init__(self, serotype: str, known_disruptions: str, disrupted_genes: str):
        self.serotype = serotype
        self.known_disruptions = known_disruptions
        self.disrupted_genes = disrupted_genes

    def read_csv(self, csv_path: str) -> pd.DataFrame:
        df = pd.read_csv(csv_path)
        return df

    def remove_known_disruptions(
        self, disrupted_df: pd.DataFrame, known_disruptions_df: pd.DataFrame
    ) -> pd.DataFrame:
        # check if serotype has a known disruption
        if self.serotype in known_disruptions_df["serotype"].values:
            # remove known disruptions from the dataframe
            genes_to_remove = known_disruptions_df.loc[
                known_disruptions_df["serotype"] == self.serotype, "gene"
            ].tolist()
            disrupted_df = disrupted_df[~disrupted_df["gene"].isin(genes_to_remove)]
        return disrupted_df

    def find_unique_samples(self, disrupted_df: pd.DataFrame) -> list:
        gene_sets = disrupted_df.groupby("sample")["gene"].apply(set).reset_index()
        # Remove duplicate gene sets and keep one representative sample per unique set
        unique_samples = gene_sets.drop_duplicates(subset=["gene"])["sample"].tolist()
        return unique_samples
