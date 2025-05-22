## Examples

This section outlines various usage scenarios for the Moncla Lab Illumina Pipeline, showcasing its flexibility and the features tested by each example configuration.

There is also a bash script `sra-fetch.sh` which uses `fasterq-dump`. See here for [more information about fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump).

### Example 1: `anna-pabird` - Large-Scale H5N1 Analysis with Local ZIP Reference

This example demonstrates a comprehensive analysis of numerous H5N1 avian influenza samples.

*   **Configuration (`anna-pabird/config.yml`):**
    *   `reference: "~/Documents/deep-sequencing/testing-data/mlip-h5n1-reference.zip"`: Utilizes a **local ZIP archive** for reference sequences. The tilde (`~`) indicates home directory expansion.
    *   `data_root_directory: "~/BaseSpace"`: Data is sourced from a standard BaseSpace download location.
*   **Metadata (`anna-pabird/metadata.tsv`):**
    *   **Extensive Dataset:** A large number of samples and associated sequencing experiments.
    *   **Varied Replicate Structures:**
        *   Many samples possess two analytical replicates.
        *   Some samples are defined with a single analytical replicate.
    *   **Composite Replicates:** Multiple `SequencingId` entries contribute to a single analytical replicate for several samples (e.g., for sample `be_w1`, `BE_W1_Rep3` and `BE_W_R1` both form `Replicate 1`).
    *   `SequencingId` values adhere to the `_Seq#` (or `_Rep#`, `_R#`) BaseSpace convention.
*   **Demonstrated Features:**
    *   **Local ZIP File as Reference Source:** Tests the pipeline's ability to unpack and utilize user-provided reference sequences from a `.zip` file.
    *   Processing of large, diverse sample sets.
    *   Handling of samples with **single and dual analytical replicates**.
    *   **Concatenation of multiple sequencing experiments** into a single analytical replicate.
    *   Standard BaseSpace data ingestion (`flow` command, default mode).

### Example 2: `anna-test` - Focused H5N1 Test with NCBI Reference

This example illustrates a smaller, targeted analysis using NCBI-derived references.

*   **Configuration (`anna-test/config.yml`):**
    *   `reference: "h5n1"`: Employs a **named reference** (`h5n1`) defined in `references.tsv`, triggering automated fetching from GenBank.
    *   `data_root_directory: "~/BaseSpace"`: Standard BaseSpace data source.
*   **Metadata (`anna-test/metadata.tsv`):**
    *   **Targeted Dataset:** Fewer samples.
    *   **Mixed Replicate Definitions:**
        *   Sample `g_com1`: Two replicates; `Replicate 2` is a composite of `G_com1_Rep2` and `G_Com1_Rep4`.
        *   Sample `kc_com1`: Two replicates, each from a distinct sequencing experiment.
        *   Sample `rf_w1`: A single analytical replicate.
*   **Demonstrated Features:**
    *   **Named Reference from `references.tsv`:** Fetches reference data directly from NCBI.
    *   Flexible handling of replicate definitions:
        *   Dual replicates with one replicate formed by **concatenating multiple sequencing experiments**.
        *   Dual replicates, each from a single sequencing experiment.
        *   **Single replicate** samples.
    *   Standard BaseSpace data ingestion.

### Example 3: `cambodia` - SRA Data Analysis with Custom Data Root

This example processes H5N1 samples sourced from the Sequence Read Archive (SRA).

*   **Configuration (`cambodia/config.yml`):**
    *   `reference: "h5n1"`: Uses a named reference from `references.tsv`.
    *   `data_root_directory: "~/Documents/deep-sequencing/testing-data/cambodia/data"`: Specifies a **custom local directory** for input FASTQ files.
*   **Metadata (`cambodia/metadata.tsv`):**
    *   **SRA Accessions:** `SequencingId` entries are SRA run accessions (e.g., `SRR9213348`).
    *   **Single Replicate per Sample:** Each SRA accession is treated as an independent sample with one analytical replicate.
*   **Demonstrated Features:**
    *   Named reference usage from `references.tsv`.
    *   Configuration of a **custom `data_root_directory`**.
    *   Ingestion of data using **SRA accessions as `SequencingId`s** (requires `python mlip/dataflow.py flow --sra-mode`).
    *   Processing of multiple samples, each with a **single analytical replicate**.

### Example 4: `stephen-canine` - Canine H3N2 Analysis

This example focuses on canine H3N2 influenza samples with a consistent dual-replicate design.

*   **Configuration (`stephen-canine/config.yml`):**
    *   `reference: h3n2`: Utilizes a **different named reference** (`h3n2`) from `references.tsv`.
    *   `data_root_directory: "~/BaseSpace"`: Standard BaseSpace data source.
*   **Metadata (`stephen-canine/metadata.tsv`):**
    *   **Uniform Replicate Structure:** All samples (`c1` to `c21`) are configured with two analytical replicates.
    *   Each replicate originates from a unique `SequencingId` (e.g., `C1` for `Replicate 1` of `c1`, `C1_seq2` for `Replicate 2`).
*   **Demonstrated Features:**
    *   Application of a **distinct named reference** for a different virus type.
    *   Consistent processing of multiple samples, all featuring **two analytical replicates**, with each replicate derived from one sequencing experiment.
    *   Standard BaseSpace data ingestion.

### Example 5: `stephen-sarscov2` - SARS-CoV-2 SRA Analysis with Modified Remapping

This example analyzes SARS-CoV-2 samples from SRA with a reduced number of remapping iterations.

*   **Configuration (`stephen-sarscov2/config.yml`):**
    *   `reference: "sarscov2"`: Uses a **SARS-CoV-2 named reference**.
    *   `data_root_directory: "~/Documents/deep-sequencing/testing-data/stephen-sarscov2/data"`: Custom data root directory.
    *   `number_of_remappings: 1`: **Remapping iterations are set to 1** (default is typically 2).
*   **Metadata (`stephen-sarscov2/metadata.tsv`):**
    *   **SRA Accessions:** `SequencingId`s are SRA run accessions.
    *   **Single Replicate per Sample:** All samples (`s1` to `s5`) are treated as single-replicate analyses.
*   **Demonstrated Features:**
    *   Use of a `sarscov2` named reference.
    *   Custom `data_root_directory`.
    *   Ingestion of SRA data (`--sra-mode` for `flow`).
    *   Processing of multiple single-replicate samples.
    *   **Customization of `number_of_remappings`**, altering the iterative refinement depth.

### Example 6: `stephen-test` - Complex H5N1 Test with Mixed Replicate Structures

This example serves as a rigorous test of the pipeline's ability to handle diverse and complex replicate definitions within a single run.

*   **Configuration (`stephen-test/config.yml`):**
    *   `reference: "h5n1"`: Named reference from `references.tsv`.
    *   `data_root_directory: "~/BaseSpace"`: Standard BaseSpace data source.
*   **Metadata (`stephen-test/metadata.tsv`):**
    *   **Highly Varied Replicate Configurations:**
        *   Samples `bv_w4`, `cb_com2`: Two replicates, each from one sequencing experiment.
        *   Samples `bv_w7`, `rf_w1`: Single replicate samples.
        *   Sample `kc_com4`: Two replicates; `Replicate 1` is a composite of `KC_com4_Rep1` and `KC_Com4_Rep3`.
        *   Sample `g_com1`: Two replicates; both `Replicate 1` (from `G_com1_Rep1`, `G_Com1_Rep3`) and `Replicate 2` (from `G_com1_Rep2`, `G_Com1_Rep4`) are composite.
*   **Demonstrated Features:**
    *   Named reference usage.
    *   Standard BaseSpace data ingestion.
    *   Robust handling of a **complex mix of replicate structures** in a single batch:
        *   Standard dual replicates.
        *   Standard single replicates.
        *   Dual replicates where one or both replicates are **composites of multiple sequencing experiments**.