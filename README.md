# Network Comparison Using Matrix-Based Ortholog Mapping

## Introduction
This code implements a matrix-based approach to compare gene coexpression networks between two species while accounting for orthologous relationships. It is particularly useful when dealing with species that have complex orthology relationships (one-to-many or many-to-one) due to gene duplication events. The method efficiently identifies conserved edges between networks by leveraging sparse matrix operations.

## Core Concept: Matrix-Based Edge Mapping
The central idea of this approach is to use matrix multiplication to efficiently map edges from one network to another through orthologous relationships. This is accomplished through three key matrices:

1. **Ortholog Matrix**: Maps genes from species A to their orthologs in species B
2. **Network Adjacency Matrix**: Represents gene relationships within each species
3. **Mapped Edges Matrix**: Result showing which edges from species A are conserved in species B

### Example to Illustrate the Process
Let's walk through a simple example with:
- 3 sorghum genes: Sb1, Sb2, Sb3
- 4 sugarcane genes: Sc1A, Sc1B (orthologs of Sb1), Sc2 (ortholog of Sb2), Sc3 (ortholog of Sb3)

#### Step 1: Ortholog Matrix (3 sorghum genes × 4 sugarcane genes)
```
     Sc1A  Sc1B  Sc2  Sc3
Sb1    1     1    0    0    # Sb1 has two orthologs: Sc1A and Sc1B
Sb2    0     0    1    0    # Sb2 has one ortholog: Sc2
Sb3    0     0    0    1    # Sb3 has one ortholog: Sc3
```

#### Step 2: Network Adjacency Matrices

Sorghum Network Matrix (3×3):
```
     Sb1  Sb2  Sb3
Sb1   0    1    0    # Sb1 connects to Sb2
Sb2   1    0    1    # Sb2 connects to Sb1 and Sb3
Sb3   0    1    0    # Sb3 connects to Sb2
```

Sugarcane Network Matrix (4×4):
```
     Sc1A  Sc1B  Sc2  Sc3
Sc1A   0     0    1    0    # Sc1A connects to Sc2
Sc1B   0     0    1    1    # Sc1B connects to Sc2 and Sc3
Sc2    1     1    0    0    # Sc2 connects to Sc1A and Sc1B
Sc3    0     1    0    0    # Sc3 connects to Sc1B
```

#### Step 3: Matrix Multiplication Result
The operation `(ortholog_matrix %*% sugarcane_matrix %*% t(ortholog_matrix))` yields:
```
     Sb1  Sb2  Sb3
Sb1   2    1    0    # Sb1 connects to Sb2 and Sb3
Sb2   1    0    0    # Sb2 connects to Sb1
Sb3   0    0    0    # Sb3 connects to Sb1
```

Comparing the original sorghum network with the mapped edges matrix, we can see which edges are conserved:
- Sb1-Sb2: Conserved (exists in both matrices)
- Sb2-Sb3: Not conserved (exists in sorghum but not in mapped edges)
- Sb1-Sb3: Does not  exist in sorghum network, neither sugarcane network

A `TRUE` (1) in the mapped edges matrix indicates that the corresponding sorghum genes have at least one pair of their orthologs connected in the sugarcane network. This demonstrates how the matrix multiplication approach can identify both conserved edges and potential new connections through ortholog relationships.

## Input Files

### 1. Orthogroups File
- **Format**: Tab-separated file generated by OrthoFinder
- **Content**: Contains ortholog relationships between species
- **Example path**: `Results_Aug30/Orthogroups/Orthogroups.tsv`

### 2. Network Edge Files
Two network files are required:
- **Format**: Tab-separated files with columns:
  - gene1: First gene in the relationship
  - gene2: Second gene in the relationship
  - weight: Edge weight (coexpression strength)
- **Example paths**:
  - Sorghum: `networks/abs_triplets_p60_sorghum_fancy.tsv`
  - Sugarcane: `networks/abs_triplets_p60_cane_fancy.tsv`

## Usage
```R
# Read and process ortholog relationships
orthologues <- read_orthogroups("path/to/Orthogroups.tsv")

# Read network files
adjacency_sorghum <- read_tsv("path/to/sorghum_network.tsv",
                             col_names = c("gene1", "gene2", "weight"))
adjacency_sugarcane <- read_tsv("path/to/sugarcane_network.tsv",
                               col_names = c("gene1", "gene2", "weight"))

# Process networks and calculate Jaccard index
results <- process_networks_matrix(adjacency_sorghum, 
                                 adjacency_sugarcane, 
                                 ortholog_pairs)
jaccard_index <- calculate_jaccard(results)
```

## Output
The code returns a data frame of conserved edges and calculates a Jaccard index measuring network similarity. The Jaccard index ranges from 0 (no similarity) to 1 (identical networks), providing a quantitative measure of network conservation between species.
