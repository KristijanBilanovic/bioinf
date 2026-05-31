#  Discovering the Gene Variants - Bioinformatics 1 project

Computational pipeline for discovering gene variants from FASTQ sequencing reads.  
Implements two approaches: **Standard (MSA-based)** and **Minimizers (fast k-mer based)**.

> Authors: Kristijan Bilanović & Luka Nikolić  
> Faculty of Electrical Engineering and Computing, University of Zagreb

---

## Requirements

- C++17 compiler (gcc ≥ 8.0 or clang ≥ 5.0)
- CMake ≥ 3.10
- Linux or macOS (Windows via WSL)

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/KristijanBilanovic/bioinf
cd bioinf
```

### 2. Set up SPOA (Multiple Sequence Alignment library)

```bash
mkdir -p vendor
cd vendor
git clone https://github.com/rvaser/spoa && cd spoa
cmake -B build -DCMAKE_BUILD_TYPE=Release
make -C build
cd ../..
```

### 3. Verify BioParser (bundled with SPOA)

```bash
ls -d ./vendor/spoa/build/_deps/bioparser-build/
ls -d ./vendor/spoa/build/_deps/bioparser-src/
ls -d ./vendor/spoa/build/_deps/bioparser-subbuild/
```

All three directories should be listed. If any are missing, redo step 2.

---

## Preparing Input Data

Place your files in the `data/` directory:

```bash
mkdir -p data/fastq/
```

Expected structure:

```
data/
├── J29B_expected.fasta      # ground truth reference (for validation)
├── J30B_expected.fasta
└── fastq/
    ├── J1.fastq
    ├── J2.fastq
    └── ...
```

- FASTQ sequencing reads go in `data/fastq/`
- Ground truth FASTA references go directly in `data/`

---

## Build

```bash
mkdir -p build
cd build
cmake ..
make
ls -la jelen_analiza   # verify the executable was created
cd ..
```

---

## Run

```bash
cd build
./jelen_analiza
```


---

## Interactive CLI

### Step 1 - File selection

```
========== File Selection ==========
1) Analyze a single file
2) Analyze all files in ../data/fastq/
Enter choice (1 or 2):
```

- **Option 1**: Analyze one FASTQ file. You'll be prompted for the filename and the ground truth reference to validate against (e.g. `J29B` or `J30B`).
- **Option 2**: Automatically processes all files starting with `J` in `data/fastq/` and runs cross-sample validation across all samples.

### Step 2 - Pipeline selection

```
========== Pipeline Selection ==========
1) Standard (parse -> filter -> msa -> cluster -> consensus)
2) With Minimizers (parse -> filter -> cluster -> consensus)
Enter choice (1 or 2):
```

| | Standard Pipeline | Minimizers Pipeline |
|---|---|---|
| **Method** | Multiple Sequence Alignment (MSA) + Hamming distance clustering | k-mer minimizers + LIS-based greedy clustering |
| **Speed** | Slower | **10–50× faster** |
| **Memory** | ~100 MB | ~10 MB |
| **Best for** | Small datasets, validation | Large datasets, batch/exploratory analysis |

---

## Output Files

After execution, results are written to:

| Output | Path |
|--------|------|
| Consensus sequences (FASTA) | `data/clusters/<SAMPLE_NAME>_CLUSTERS.fasta` |
| Cluster membership records | `data/cluster_members/<SAMPLE_NAME>_CLUSTER_MEMBERS.txt` |
| Cross-sample validation report | `data/clusters/cross_sample_analysis.txt` |

---

## Project Structure

```
bioinf/
├── CMakeLists.txt
├── main.cpp
├── include/
│   └── sequence_analyzer.hpp
├── src/
│   └── sequence_analyzer.cpp
├── data/
│   └── fastq/
└── vendor/
    └── spoa/
```

---

## Key Parameters

```cpp
// Minimizers pipeline
K_MER_LENGTH      = 11     // k-mer size
WINDOW_SIZE       = 5      // minimizer window
DISTANCE_THRESHOLD = 0.33  // clustering threshold (≤ 0.33 → same cluster)

// Standard pipeline (MSA + Hamming)
ALIGNMENT_MATCH    =  3
ALIGNMENT_MISMATCH = -5
ALIGNMENT_GAP      = -3
HAMMING_CLUSTER_K  = 12    // max Hamming distance to merge aligned reads

// Cross-sample validation
LARGE_CLUSTER_THRESHOLD = 100   // reads needed to count as a "large" cluster
HAMMING_THRESHOLD       = 5     // max distance for a cross-sample match
```

---

## References

1. Križanović, K. *Bioinformatika 1.* FER, University of Zagreb.
2. Vaser, R. et al. *SPOA: Fast and Accurate Whole-Genome Multiple Sequence Alignment.* — [github.com/rvaser/spoa](https://github.com/rvaser/spoa)
3. Kosier, S. *Pronalaženje varijanti gena iz podataka dobivenih sekvenciranjem.* FER, Zagreb, 2019.
