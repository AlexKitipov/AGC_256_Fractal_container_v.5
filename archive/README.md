# AGC_256_Fractal_container_v.5
“AGC‑128/256 Notepad is an experimental framework for DNA‑inspired data encoding, Unicode/ASCII transformation, FASTA metadata handling, genetic checksums, and conceptual AGC‑256 fractal logic. Includes a modular architecture and optional GUI for exploration and testing.”

"""## AGC-256: The Five Laws of Fractal DNA Encoding

This document outlines the foundational principles, or 'laws,' that govern the Adaptive Genetic Code 256 (AGC-256) system. These laws describe how bit streams are structured, how nucleotides connect, the fractal nature of its containers, and how larger 'cubes' are derived from smaller components. The Python functions mentioned implement these conceptual laws.

---

### I. Law of Sliding Windows: Bit Stream Analysis

This law defines how overlapping windows of varying bit lengths are analyzed within a given bit stream. It's fundamental for understanding local structural patterns.

**Concept**: For a bit stream of length `n`, an overlapping window of size `k` (where `k` is the window length) shifts one bit at a time. The number of such windows is directly related to the stream length and window size.

**Formulas**:
*   **2-bit windows (Nucleotide)**: `n - 2 + 1`
*   **3-bit windows (Triple Motif)**: `n - 3 + 1`
*   **4-bit windows (Cube)**: `n - 4 + 1`

**Implemented Functions**:
*   `calculate_2_bit_windows(bit_stream_length)`
*   `calculate_3_bit_windows(bit_stream_length)`
*   `calculate_4_bit_windows(bit_stream_length)`

**Example**:
For a bit stream length of 10:
*   Number of 2-bit windows: 9
*   Number of 3-bit windows: 8
*   Number of 4-bit windows: 7

---

### II. Law of Nucleotide Connections

This law describes the combinatorial relationships between individual nucleotides within a sequence, crucial for understanding sequence complexity and potential interactions.

**Concept**: Given a set of `n` nucleotides, one can calculate the number of unique ordered pairs and the total number of connections when considering different interaction 'modes'.

**Formulas**:
*   **Number of ordered nucleotide pairs**: `(n * (n - 1)) / 2`
*   **Number of connections with 2 modes**: `n * (n - 1)`

**Implemented Functions**:
*   `calculate_ordered_nucleotide_pairs(num_nucleotides)`
*   `calculate_connections_with_2_modes(num_nucleotides)`

**Example**:
For 4 nucleotides:
*   Number of ordered nucleotide pairs: 6.0
*   Number of connections with 2 modes: 12

---

### III. Law of Nested Brackets (Fractal Container)

This law introduces the conceptual `FractalCube` as the fundamental unit of information in AGC-256, characterized by a nested, self-similar structure of information layers.

**Concept**: An AGC-256 'Cube' is a fractal container comprising a central 'core' (a 2-bit motif or nucleotide), an 'internal context' (flanking 1-bit motifs directly adjacent to the core), and an 'external context' (outermost 1-bit motifs providing broader context to the internal layer). This nested hierarchy allows for multi-scale information embedding.

**Implemented Class**:
*   `FractalCube`
    *   `core`: The central 2-bit motif (e.g., 'A', 'T', 'G', 'C').
    *   `internal_context`: A tuple of two 1-bit motifs flanking the core (e.g., `('0', '1')`).
    *   `external_context`: A tuple of two 1-bit motifs providing outer context (e.g., `('0', '1')`).

**Example**:
Creating an instance of `FractalCube`:
```
FractalCube(
  core='AG',
  internal_context=('0', '1'),
  external_context=('0', '0'))
```

---

### IV. Law of Fractality

This law quantifies the self-similar patterns or 'motives' that can be found within a bit stream, highlighting the fractal nature of AGC-256 information.

**Concept**: The number of distinct motives (patterns) of a certain size `k` that can be generated or identified within a bit stream of length `n` follows a specific fractal relationship. This indicates the rich structural information embedded at various scales.

**Formula**:
*   Number of motives `F_k` for a motive size `k` in a bit stream of length `n`: `(n - k + 1) * (2 ** k)`

**Implemented Function**:
*   `calculate_fractal_motives(bit_stream_length, motive_size)`

**Example**:
For a bit stream length of 10:
*   Number of motives of size 2: 36
*   Number of motives of size 3: 64

---

### V. Law of AGC-128 → AGC-256 Derivation

This law describes the mechanistic process by which a 4-bit AGC-256 'Cube' is constructed from a 2-bit AGC-128 core and contextual single bits, representing the fundamental building block transition from AGC-128 to AGC-256.

**Concept**: An AGC-256 cube is a higher-order structure that integrates a 2-bit AGC-128 core nucleotide with a 1-bit context from its left and a 1-bit context from its right. These elements combine to form a complete 4-bit unit, where the internal and external contexts are derived from these flanking bits.

**Mechanism**: The 2-bit core nucleotide is converted to its binary representation. This 2-bit binary string is then sandwiched between the left and right 1-bit context bits to form a 4-bit binary string. This 4-bit string defines the structure of the `FractalCube`, where the internal and external contexts are directly represented by the original context bits.

**Implemented Function**:
*   `simulate_agc256_cube_derivation(core_nucleotide, left_context_bit, right_context_bit)`

**Example**:
Simulating AGC-256 Cube derivation from `core_nucleotide='A'`, `left_context_bit='0'`, `right_context_bit='1'`:
```
FractalCube(
  core='A',
  internal_context=('0', '1'),
  external_context=('0', '1'))
```

# AGC_256_Fractal_container_v.5 - Adaptive Genetic Code (Unified & Fractal Edition)

### Official README (Unified & Fractal Edition)

## Authors
- **Aleksandar Kitipov**  
  Emails: aeksandar.kitipov@gmail.com / aeksandar.kitipov@outlook.com  
- **Copilot**  
  Co‑author, technical collaborator, documentation support
  - **Gemini 2.5 Flash**  
  Co‑author, technical collaborator, documentation support

---

## 1. Overview
**AGC_256_Fractal_container_v.5** is the unified and enhanced version of the Adaptive Genetic Code system, significantly expanding its capabilities from text encoding to universal binary data and conceptualizing a fractal architecture (AGC-256). It combines the original **v.1 (ASCII)** capabilities with the **v.2 (Unicode)** extension, adds robust **binary file handling**, and integrates a **self-recovery mechanism**. Furthermore, it lays the theoretical groundwork for **AGC-256** by implementing functions for its foundational mathematical and structural 'laws'.

This version features an interactive Graphical User Interface (GUI) built with `tkinter`, allowing users to seamlessly switch between v1 and v2 encoding/decoding, handle binary files, manage genetic checksums, and explore self-recovery options. Like its predecessors, it requires **no external libraries** for its core encoding/decoding logic, maintaining a tiny footprint while ensuring high data integrity through its self-checking genetic structure.

---

## 2. What the Program Does

AGC_256_Fractal_container_v.5 provides complete reversible transformations:

### v.1 (ASCII) Transformation:
```
Text → ASCII (8 bits) → 4 (2-bit genes) → A/T/G/C DNA Sequence
```
And back:
```
DNA Sequence → 4 (2-bit genes) → 8-bit ASCII → Text
```

### v.2 (Unicode) Transformation:
```
Unicode Text → UTF-8 Bytes → Length Gene + Genetic Bytes → A/T/G/C DNA Sequence
```
And back:
```
DNA Sequence → Genetic Bytes + Length Gene → UTF-8 Bytes → Unicode Text
```

### Binary File Transformation:
```
Binary File → Raw Bytes → 4 (2-bit genes per byte) → A/T/G/C DNA Sequence (with metadata in FASTA header)
```
And back:
```
DNA Sequence (with metadata) → 4 (2-bit genes per byte) → Raw Bytes → Binary File
```

This system precisely preserves:
- **v.1:** Letters, numbers, punctuation, whitespace, and ASCII extended symbols.
- **v.2:** Any Unicode character (including ASCII, Cyrillic, CJK, Emojis, Symbols), preserving structured blocks and FASTA-formatted sequences.
- **Binary:** Any arbitrary binary data (images, executables, compressed files, etc.), along with its original filename, extension, and size.

**If you encode data and decode it again using the correct version, the output will match the original exactly, character-for-character or byte-for-byte.**

---

## 3. Key Features

### 3.1. Unified Encoding/Decoding
Supports `v1 (ASCII)`, `v2 (Unicode)`, and **Binary File** encoding/decoding for a universal data container approach.

### 3.2. Intuitive Graphical User Interface (GUI)
Built with `tkinter` for ease of use, featuring:
- **Version Selection:** Radio buttons to choose between `v1 (ASCII)` and `v2 (Unicode)` modes.
- **File Operations:** Open, Save, Save As, and New file functionalities.
- **Edit Menu:** Standard text editor actions like Undo, Redo, Cut, Copy, Paste, Delete, and Select All.
- **Context Menu:** Right-click menu for quick editing actions.
- **FASTA Management:** Encode text or binary files to FASTA and Load/Decode text or binary FASTA files.
- **Checksum Integration:** Options to add, verify, and consider genetic checksums during encoding/decoding.
- **Visualization Placeholder:** Provides basic sequence information.
- **Binary Representation Viewer**: Displays the raw binary string of the current nucleotide sequence.

### 3.3. Full Reversibility
Every piece of data, whether ASCII, Unicode, or binary, is transformed reversibly, ensuring zero data loss upon decoding.

### 3.4. Self-Checking Genetic Structure
AGC-128 maintains its three core biological-style integrity rules:
- **Sum-2 Rule**: Each 2-bit gene has a total bit-sum of 2. Any bit flip breaks the rule and becomes detectable.
- **No-Triple Rule**: The sequence can never contain `111` or `000`. If such a pattern appears, the data is invalid.
- **Deterministic-Next-Bit Rule**: Predictable bit sequences (`11` → `0`, `00` → `1`). This allows partial reconstruction of missing or damaged data.

### 3.5. FASTA Compatibility
The DNA output can be saved as a `.fasta` file, complete with embedded metadata for binary files, making it suitable for digital archiving, DNA-like storage experiments, and bioinformatics-style workflows.

### 3.6. Metadata Handling for Binary Files
Integrates functions to extract (filename, extension, size, type=BINARY) and embed metadata into the FASTA header during encoding, and to parse and utilize this metadata for accurate file reconstruction during decoding.

### 3.7. Self-Recovery Capabilities
Implements functions for:
- **Error Detection**: Identifies violations of Sum-2, No-Triple, and Deterministic-Next-Bit rules.
- **Deterministic Reconstruction**: Corrects small, unambiguous errors based on rules.
- **Conceptual Variant Generation/Selection**: Provides a framework for generating and selecting alternative valid sequences for larger corrupted segments using checksums.
- **Recovery Confidence**: Calculates a score (0-100) indicating the reliability of the recovery process.

### 3.8. AGC-256 Conceptual Framework
Incorporates functions and classes demonstrating the theoretical 'laws' that govern a potential future AGC-256 system, focusing on fractal bit stream analysis and container structures.

---

## 4. Genetic Alphabet
AGC-128 uses four genetic symbols mapped from 2-bit pairs:

```
11 → G  
00 → C  
10 → A  
01 → T
```

---

## 5. AGC-128 v2 (Unicode) Core Principles in Detail

### 5.1. UTF-8 as Foundation
Unicode characters are first converted to their UTF-8 byte representation (1 to 4 bytes).

### 5.2. Length Prefix Gene
Each encoded Unicode character begins with a single-nucleotide `Length Gene` that indicates the number of UTF-8 bytes that follow for that character:

| UTF-8 Length | Number of Bytes | 2-bit Marker | Length Gene |
|--------------|-----------------|--------------|-------------|
| 1 byte       | ASCII           | 00           | C           |
| 2 bytes      | Cyrillic        | 01           | T           |
| 3 bytes      | Multi-byte      | 10           | A           |
| 4 bytes      | Emojis          | 11           | G           |

### 5.3. Byte Encoding
Each individual UTF-8 byte (0-255) is encoded into four 2-bit nucleotide genes, consistent with AGC-128 v1's 8-bit to 4-nucleotide conversion.

Thus, a Unicode character's genetic sequence is: `[Length Gene] + [4 genes per byte]`.
- **1-byte UTF-8 (ASCII)** → `C` + 4 genes = 5 nucleotides
- **2-bytes UTF-8 (e.g., Cyrillic)** → `T` + 8 genes = 9 nucleotides
- **3-bytes UTF-8 (e.g., Chinese)** → `A` + 12 genes = 13 nucleotides
- **4-bytes UTF-8 (e.g., Emojis)** → `G` + 16 genes = 17 nucleotides

---

## 6. Binary File Handling

For binary files, the content is read directly as raw bytes. Each byte is converted into 4 nucleotides. Metadata (original filename, extension, size) is extracted and serialized into the FASTA header. During decoding, this metadata is parsed to reconstruct the original file, ensuring that all binary data (e.g., images, audio, executables) can be perfectly restored.

---

## 7. Genetic Checksum

An optional 2-nucleotide genetic checksum can be appended to the entire sequence to verify data integrity. It calculates the sum of all 2-bit nucleotide values, modulo 16, and encodes this 4-bit result into two nucleotides. The GUI provides explicit options to add and verify this checksum, ensuring flexibility and data validation. For binary files, if a checksum is expected, it will be removed before binary decoding, even if found to be invalid, to allow file reconstruction (with appropriate warnings about data integrity).

---

## 8. Self-Recovery System

The self-recovering AGC-128 system is capable of detecting and, in some cases, automatically correcting errors in nucleotide sequences based on predefined rules:

-   **'Sum-2 Rule'**: Ensures all nucleotides are valid ('A', 'T', 'G', 'C').
-   **'No-Triple Rule'**: Checks for absence of '000' or '111' binary patterns.
-   **'Deterministic-Next-Bit Rule'**: Verifies specific bit sequences ('11' followed by '0', and '00' followed by '1').

It can perform deterministic single-bit corrections for 'Deterministic-Next-Bit Rule' violations. A framework for generating and evaluating candidate variants for more complex errors using a checksum-based selection process is also conceptualized. The system provides a confidence score (0-100) indicating the success and reliability of the recovery process.

---

## 9. AGC-256: The Five Laws of Fractal DNA Encoding

This section outlines the foundational principles, or 'laws,' that govern the Adaptive Genetic Code 256 (AGC-256) system. These laws describe how bit streams are structured, how nucleotides connect, the fractal nature of its containers, and how larger 'cubes' are derived from smaller components. The Python functions mentioned implement these conceptual laws.

### I. Law of Sliding Windows: Bit Stream Analysis

This law defines how overlapping windows of varying bit lengths are analyzed within a given bit stream. It's fundamental for understanding local structural patterns.

**Concept**: For a bit stream of length `n`, an overlapping window of size `k` (where `k` is the window length) shifts one bit at a time. The number of such windows is directly related to the stream length and window size.

**Formulas**:
*   **2-bit windows (Nucleotide)**: `n - 2 + 1`
*   **3-bit windows (Triple Motif)**: `n - 3 + 1`
*   **4-bit windows (Cube)**: `n - 4 + 1`

**Implemented Functions**:
*   `calculate_2_bit_windows(bit_stream_length)`
*   `calculate_3_bit_windows(bit_stream_length)`
*   `calculate_4_bit_windows(bit_stream_length)`

### II. Law of Nucleotide Connections

This law describes the combinatorial relationships between individual nucleotides within a sequence, crucial for understanding sequence complexity and potential interactions.

**Concept**: Given a set of `n` nucleotides, one can calculate the number of unique ordered pairs and the total number of connections when considering different interaction 'modes'.

**Formulas**:
*   **Number of ordered nucleotide pairs**: `(n * (n - 1)) / 2`
*   **Number of connections with 2 modes**: `n * (n - 1)`

**Implemented Functions**:
*   `calculate_ordered_nucleotide_pairs(num_nucleotides)`
*   `calculate_connections_with_2_modes(num_nucleotides)`

### III. Law of Nested Brackets (Fractal Container)

This law introduces the conceptual `FractalCube` as the fundamental unit of information in AGC-256, characterized by a nested, self-similar structure of information layers.

**Concept**: An AGC-256 'Cube' is a fractal container comprising a central 'core' (a 2-bit motif or nucleotide), an 'internal context' (flanking 1-bit motifs directly adjacent to the core), and an 'external context' (outermost 1-bit motifs providing broader context to the internal layer). This nested hierarchy allows for multi-scale information embedding.

**Implemented Class**:
*   `FractalCube`
    *   `core`: The central 2-bit motif (e.g., 'A', 'T', 'G', 'C').
    *   `internal_context`: A tuple of two 1-bit motifs flanking the core (e.g., `('0', '1')`).
    *   `external_context`: A tuple of two 1-bit motifs providing outer context (e.g., `('0', '1')`).

### IV. Law of Fractality

This law quantifies the self-similar patterns or 'motives' that can be found within a bit stream, highlighting the fractal nature of AGC-256 information.

**Concept**: The number of distinct motives (patterns) of a certain size `k` that can be generated or identified within a bit stream of length `n` follows a specific fractal relationship. This indicates the rich structural information embedded at various scales.

**Formula**:
*   Number of motives `F_k` for a motive size `k` in a bit stream of length `n`: `(n - k + 1) * (2 ** k)`

**Implemented Function**:
*   `calculate_fractal_motives(bit_stream_length, motive_size)`

### V. Law of AGC-128 → AGC-256 Derivation

This law describes the mechanistic process by which a 4-bit AGC-256 'Cube' is constructed from a 2-bit AGC-128 core and contextual single bits, representing the fundamental building block transition from AGC-128 to AGC-256.

**Concept**: An AGC-256 cube is a higher-order structure that integrates a 2-bit AGC-128 core nucleotide with a 1-bit context from its left and a 1-bit context from its right. These elements combine to form a complete 4-bit unit, where the internal and external contexts are derived from these flanking bits.

**Mechanism**: The 2-bit core nucleotide is converted to its binary representation. This 2-bit binary string is then sandwiched between the left and right 1-bit context bits to form a 4-bit binary string. This 4-bit string defines the structure of the `FractalCube`, where the internal and external contexts are directly represented by the original context bits.

**Implemented Function**:
*   `simulate_agc256_cube_derivation(core_nucleotide, left_context_bit, right_context_bit)`

---

## 10. Usage

### Local GUI Execution:
To use the GUI, save the provided consolidated Python code as a `.py` file (e.g., `agc_notepad.py`). Then, run it from your local machine's terminal using a Python interpreter that supports Tkinter (typically pre-installed with Python on most desktop OSes):

```bash
python agc_notepad.py
```

**GUI Features:**
1.  **Select Encoding/Decoding Version:** Use the radio buttons (`v1 (ASCII)` or `v2 (Unicode)`).
2.  **Type or Load Text:** Enter text directly or use `File > Open`.
3.  **Encode Text to FASTA:** Go to `Encode > Encode Text to AGC-128 FASTA`. You'll be prompted for a FASTA header and whether to add a checksum.
4.  **Encode Binary File to FASTA:** Go to `Encode > Binary > Encode Binary File to AGC-128 FASTA`. Select a binary file; metadata will be automatically embedded.
5.  **Load and Decode Text FASTA:** Go to `Decode > Load and Decode Text FASTA`. The system will prompt if a checksum is expected and verify it if present.
6.  **Load and Decode Binary FASTA:** Go to `Decode > Binary > Decode AGC-128 to Binary File`. The system will parse metadata and save the reconstructed binary file to its original format.
7.  **Tools:** `Verify Checksum` (for currently loaded sequence), `Visualize Sequence` (placeholder info), `View Binary Representation`, and `Attempt Self-Recovery`.

### Command-Line Usage (Colab/Non-GUI Environment):
When running in environments without a GUI (like Google Colab), the script automatically detects this and will instead print console demonstrations of its core encoding/decoding, checksum, self-recovery, and AGC-256 conceptual functions.

---

## 11. Project Status
- **AGC_256_Fractal_container_v.5** — Stable unified core with GUI, binary support, self-recovery, and AGC-256 conceptual framework.

---

## 12. Notes
This README represents the comprehensive documentation for the unified AGC-128 Notepad, incorporating all features and improvements developed through collaborative discussions and testing. Its modular design and clear separation of concerns facilitate future expansions, such as integrating the conceptual AGC-256 functions into the GUI or enhancing the self-recovery mechanism with actual variant generation and selection.

## Final Task

### Subtask:
Summarize the capabilities of the final consolidated script and provide instructions for local execution and testing.

## Summary:

### Q&A
The final consolidated script encompasses a wide range of functionalities for the AGC-128 Notepad, including core encoding and decoding capabilities for both v1 ASCII and v2 Unicode, as well as binary file handling. It integrates robust metadata management (extraction, serialization, and parsing) for FASTA files, and implements genetic checksums for data integrity. A sophisticated self-recovery mechanism is included, featuring violation detection, deterministic reconstruction, candidate variant generation (conceptual), and a confidence scoring system. Furthermore, it incorporates conceptual functions and classes related to the AGC-256 standard, such as windowing, nucleotide connection rules, the FractalCube data structure, fractal motive calculations, and AGC-256 Cube derivation. The script also provides a Tkinter-based Graphical User Interface (GUI) for interactive use.

For local execution and testing, save the provided consolidated Python code as a `.py` file (e.g., `agc128_notepad.py`). Then, run it from your local machine's terminal using a Python interpreter that supports Tkinter (typically pre-installed with Python on most desktop OSes): `python agc128_notepad.py`. If executed in a non-GUI environment like Google Colab, the script automatically detects this and will instead print console demonstrations of its core encoding/decoding, checksum, self-recovery, and AGC-256 conceptual functions. The GUI provides menu options for file operations, encoding/decoding text or binary files to/from FASTA, verifying checksums, visualizing sequences (placeholder), viewing binary representations, and attempting self-recovery.

### Data Analysis Key Findings
*   The script successfully detects the execution environment (e.g., Google Colab), automatically switching from GUI mode to console demonstration for core functionalities.
*   **V2 Unicode Encoding/Decoding:** Demonstrated accurate encoding of complex Unicode strings (e.g., "Здравейте, свят!☺ 123") into nucleotide sequences and their subsequent lossless decoding, confirming full Unicode support.
*   **Checksum Functionality:** `add_genetic_checksum` and `verify_genetic_checksum` functions correctly applied and validated checksums for both V1 ASCII and V2 Unicode encoded sequences.
*   **AGC-256 Conceptual Functions:** All five laws of AGC-256 (windowing, nucleotide connection, FractalCube structure, fractal motives, and cube derivation) were successfully demonstrated, yielding expected numerical results and structural representations in the console output.
*   **V1 ASCII Encoding/Decoding:** Confirmed correct encoding and decoding of standard ASCII text (e.g., "Hello, Colab!") into and from nucleotide sequences.
*   **Code Structure:** The consolidated script adheres to the specified architectural requirements, with all imports at the beginning, logical function grouping, and the GUI entry point (`setup_gui()`) correctly placed within an `if __name__ == "__main__":` block.

### Insights or Next Steps
*   The modular design and clear separation of concerns within the consolidated script facilitate future expansions, such as integrating the conceptual AGC-256 functions into the GUI or enhancing the self-recovery mechanism with actual variant generation and selection.
*   Further development should focus on implementing the advanced self-recovery features, specifically the `generate_candidate_variants` and `select_best_variant_with_checksum` functions, to move beyond their current conceptual placeholder status and enable actual data repair.
"""
