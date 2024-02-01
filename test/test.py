def merge_count_files(run) -> None:
    """
    Merge count files into a single count matrix TSV file.

    Args:
        run (str): The run number.

    Returns:
        None
    """

    if run == "1":
        import os
        import pandas as pd
        import re

        # Define the directory containing the count files
        directory = "miRNA_pipe/test/"

        # Get a list of all count files in the directory
        count_files = sorted(
            [
                file
                for file in os.listdir(directory)
                if file.endswith("_miRNA_concat.txt")
            ]
        )

        # Initialize an empty DataFrame to store the merged data
        merged_data = pd.DataFrame(columns=["miRNA"])

        # Loop through each count file
        for file in count_files:
            # Read the count file into a DataFrame
            df = pd.read_csv(os.path.join(directory, file), delimiter="\t")

            # Extract the filename (excluding the file extension) to use as column header
            filename = re.sub(r"_miRNA_concat\.txt$", "", file)

            # Rename the columns, excluding the first column (miRNA)
            df.columns = ["miRNA"] + [filename for col in df.columns[1:]]

            # Merge the DataFrame with the merged_data DataFrame, using the miRNA column as the key
            merged_data = pd.merge(merged_data, df, on="miRNA", how="outer")

        # Convert all non-integer values to NaN and then replace NaN with 0
        for col in merged_data.columns[1:]:
            merged_data[col] = (
                merged_data[col]
                .apply(pd.to_numeric, errors="coerce")
                .fillna(0)
                .astype(int)
            )

        # Write the merged data to a new TSV file
        merged_data.to_csv(
            os.path.join(directory, "count_matrix.tsv"), sep="\t", index=False
        )


# Merges all the sample counts to create the count matrix.
merge_count_files(run="1")
