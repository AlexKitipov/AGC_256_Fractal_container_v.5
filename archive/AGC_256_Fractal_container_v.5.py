# AGC_256_Fractal_container_v.5
import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox, scrolledtext
import sys
import os
import copy # For deepcopy in self-recovery functions

# =========================
# GLOBAL STATE
# =========================
current_encoded_nucleotide_sequence = []

# =========================
# AGC-128 CORE TABLES
# =========================

# 00 -> C, 01 -> T, 10 -> A, 11 -> G
nuc_to_int = {
    'C': 0,
    'T': 1,
    'A': 2,
    'G': 3
}
int_to_nuc = {v: k for k, v in nuc_to_int.items()}

# For V2 Unicode
LENGTH_MAP = {
    1: 'C',  # 1 byte UTF-8 (ASCII)
    2: 'T',  # 2 bytes UTF-8 (e.g., Cyrillic)
    3: 'A',  # 3 bytes UTF-8 (other multi-byte)
    4: 'G'   # 4 bytes UTF-8 (emojis)
}
REV_LENGTH_MAP = {v: k for k, v in LENGTH_MAP.items()}

# Map 2-bit strings to nucleotides for V2 byte-level encoding
bit_to_nuc = {
    '00': 'C',
    '01': 'T',
    '10': 'A',
    '11': 'G'
}

# =========================
# ENCODING FUNCTIONS
# =========================

# V1 ASCII Encoding
def string_to_nucleotide_sequence_v1(text):
    """
    Всеки символ -> ASCII (8 бита) -> 4 двойки бита -> 4 нуклеотида.
    """
    seq = []
    for ch in text:
        ascii_val = ord(ch)
        # Extract 2-bit chunks
        b1 = (ascii_val >> 6) & 0b11  # Most significant 2 bits
        b2 = (ascii_val >> 4) & 0b11
        b3 = (ascii_val >> 2) & 0b11
        b4 = ascii_val & 0b11        # Least significant 2 bits
        seq.extend([
            int_to_nuc[b1],
            int_to_nuc[b2],
            int_to_nuc[b3],
            int_to_nuc[b4]
        ])
    return seq

# V2 Unicode Helper Functions (byte-level)
def byte_to_tagc_v2(byte):
    """
    Converts a single byte (0-255) into its corresponding 4 TAGC nucleotides.
    """
    bits = f"{byte:08b}"
    tagc_nucleotides = []
    for i in range(0, 8, 2):
        two_bit_chunk = bits[i:i+2]
        tagc_nucleotides.append(bit_to_nuc[two_bit_chunk])
    return tagc_nucleotides

# V2 Unicode Encoding
def encode_unicode_char_to_tagc(unicode_char):
    """
    Converts a single Unicode character into a TAGC nucleotide sequence,
    prefixed with a Length Gene.
    """
    utf8_bytes = unicode_char.encode('utf-8')
    num_bytes = len(utf8_bytes)
    encoded_sequence = []

    if num_bytes not in LENGTH_MAP:
        raise ValueError(f"Unsupported UTF-8 byte length: {num_bytes} for character '{unicode_char}'")

    length_gene = LENGTH_MAP[num_bytes]
    encoded_sequence.append(length_gene)

    for byte_val in utf8_bytes:
        tagc_nucleotides = byte_to_tagc_v2(byte_val)
        encoded_sequence.extend(tagc_nucleotides)

    return encoded_sequence

def encode_string_to_unicode_tagc_sequence(input_string):
    """
    Encodes an entire string into a Unicode TAGC nucleotide sequence.
    """
    full_tagc_sequence = []
    for char in input_string:
        char_tagc = encode_unicode_char_to_tagc(char)
        full_tagc_sequence.extend(char_tagc)
    return full_tagc_sequence

# Binary Encoding
def bytes_to_nucleotide_sequence(raw_bytes):
    """
    Converts a sequence of raw bytes into its AGC-128 nucleotide representation.
    Each byte is converted into 4 nucleotides.

    Args:
        raw_bytes (bytes): A bytes object (e.g., from reading a binary file).

    Returns:
        list: A list of nucleotide characters (e.g., ['A', 'T', 'G', 'C']) representing the encoded bytes.
    """
    nucleotide_sequence = []
    for byte_val in raw_bytes:
        tagc_nucleotides = byte_to_tagc_v2(byte_val)
        nucleotide_sequence.extend(tagc_nucleotides)
    return nucleotide_sequence

# =========================
# CHECKSUM FUNCTIONS
# =========================

def calculate_genetic_checksum(nucleotide_sequence):
    """
    Calculates a genetic checksum for a given nucleotide sequence.
    The checksum is based on the sum of 2-bit integer representations
    of nucleotides, modulo 16, encoded as two nucleotides.
    """
    total_sum = 0
    for nuc in nucleotide_sequence:
        total_sum += nuc_to_int.get(nuc, 0)  # Use .get with default 0 for safety

    checksum_value = total_sum % 16  # Checksum is a value between 0 and 15 (4-bit value)

    # Convert checksum value to 4-bit binary string (e.g., 0 -> "0000", 15 -> "1111")
    checksum_binary = f"{checksum_value:04b}"

    # Convert 4-bit binary string to two nucleotides using int_to_nuc
    checksum_nuc1_int = int(checksum_binary[0:2], 2)
    checksum_nuc2_int = int(checksum_binary[2:4], 2)

    checksum_nuc1 = int_to_nuc[checksum_nuc1_int]
    checksum_nuc2 = int_to_nuc[checksum_nuc2_int]

    return [checksum_nuc1, checksum_nuc2]

def add_genetic_checksum(seq):
    """
    Appends the calculated genetic checksum to a copy of the original nucleotide sequence.
    """
    checksum = calculate_genetic_checksum(seq)
    sequence_with_checksum = list(seq)  # Create a copy
    sequence_with_checksum.extend(checksum)
    return sequence_with_checksum

def verify_genetic_checksum(seq):
    """
    Verifies the genetic checksum of a sequence.
    Assumes the last two nucleotides are the checksum.
    """
    if len(seq) < 2:
        return False
    data = seq[:-2]        # The original data part
    checksum = seq[-2:]    # The provided checksum part
    expected = calculate_genetic_checksum(data)
    return checksum == expected

# =========================
# DECODING FUNCTIONS
# =========================

# V1 ASCII Decoding
def decode_nucleotide_sequence_to_string_v1(nucleotide_sequence):
    """
    4 нуклеотида -> 4x2 бита -> 8-битов ASCII.
    """
    decoded_chars = []
    for i in range(0, len(nucleotide_sequence), 4):
        chunk = nucleotide_sequence[i:i+4]
        if len(chunk) != 4:
            # Warning already handled in GUI if length mismatch
            break

        # Convert each nucleotide to its 2-bit integer representation
        b1 = nuc_to_int[chunk[0]]
        b2 = nuc_to_int[chunk[1]]
        b3 = nuc_to_int[chunk[2]]
        b4 = nuc_to_int[chunk[3]]

        # Combine the four 2-bit integers to form a single 8-bit integer
        ascii_val = (b1 << 6) | (b2 << 4) | (b3 << 2) | b4
        decoded_chars.append(chr(ascii_val))
    return "".join(decoded_chars)

# V2 Unicode Helper Functions (byte-level)
def tagc_to_byte_v2(nucleotides):
    """
    Converts 4 TAGC nucleotides back into a single byte.
    """
    if len(nucleotides) != 4:
        raise ValueError("Input must be a list of exactly 4 nucleotides.")

    binary_string = ""
    for nuc in nucleotides:
        int_value = nuc_to_int[nuc]
        binary_string += f"{int_value:02b}"

    byte_value = int(binary_string, 2)
    return byte_value

# V2 Unicode Decoding
def decode_tagc_to_unicode_char(tagc_sequence_chunk):
    """
    Decodes a chunk of TAGC nucleotides representing a single encoded Unicode character
    back into the original Unicode character.
    """
    if not tagc_sequence_chunk:
        raise ValueError("Input tagc_sequence_chunk cannot be empty.")

    length_gene = tagc_sequence_chunk[0]

    if length_gene not in REV_LENGTH_MAP:
        raise ValueError(f"Invalid Length Gene '{length_gene}' found.")
    num_bytes = REV_LENGTH_MAP[length_gene]

    expected_length = 1 + (num_bytes * 4)

    if len(tagc_sequence_chunk) != expected_length:
        raise ValueError(
            f"Mismatch in TAGC sequence chunk length. Expected {expected_length} nucleotides "
            f"but got {len(tagc_sequence_chunk)}. (Length Gene: {length_gene}, num_bytes: {num_bytes}) "
            f"Full chunk: {tagc_sequence_chunk}"
        )

    data_nucleotides = tagc_sequence_chunk[1:]
    byte_array = bytearray()

    for i in range(0, len(data_nucleotides), 4):
        nuc_chunk = data_nucleotides[i:i+4]
        decoded_byte = tagc_to_byte_v2(nuc_chunk)
        byte_array.append(decoded_byte)

    decoded_char = byte_array.decode('utf-8')
    return decoded_char

def decode_unicode_tagc_sequence_to_string(tagc_sequence):
    """
    Decodes an entire Unicode TAGC nucleotide sequence back into a string.
    """
    decoded_chars = []
    current_index = 0

    while current_index < len(tagc_sequence):
        length_gene = tagc_sequence[current_index]

        if length_gene not in REV_LENGTH_MAP:
            raise ValueError(f"Invalid Length Gene '{length_gene}' at index {current_index}.")
        num_bytes = REV_LENGTH_MAP[length_gene]

        char_chunk_length = 1 + (num_bytes * 4)

        char_tagc_chunk = tagc_sequence[current_index:current_index + char_chunk_length]

        if len(char_tagc_chunk) != char_chunk_length:
            raise ValueError(
                f"Incomplete TAGC sequence at index {current_index}. "
                f"Expected {char_chunk_length} nucleotides, but found {len(char_tagc_chunk)}."
            )

        decoded_char = decode_tagc_to_unicode_char(char_tagc_chunk)
        decoded_chars.append(decoded_char)

        current_index += char_chunk_length

    return "".join(decoded_chars)

# Binary Decoding
def nucleotide_sequence_to_bytes(nucleotide_sequence):
    """
    Converts an AGC-128 nucleotide sequence back into a raw byte sequence.
    Each 4 nucleotides are converted back into a single byte.

    Args:
        nucleotide_sequence (list): A list of nucleotide characters (e.g., ['A', 'T', 'G', 'C']).

    Returns:
        bytes: A bytes object representing the decoded data.
    """
    byte_array = bytearray()
    for i in range(0, len(nucleotide_sequence), 4):
        chunk = nucleotide_sequence[i:i+4]
        if len(chunk) != 4:
            raise ValueError(f"Incomplete nucleotide chunk for byte conversion at index {i}. Expected 4, got {len(chunk)}.")

        decoded_byte = tagc_to_byte_v2(chunk)
        byte_array.append(decoded_byte)

    return bytes(byte_array)

# =========================
# METADATA HANDLING FUNCTIONS
# =========================
def extract_file_metadata(file_path):
    """
    Extracts metadata (original filename, extension, size, type=BINARY) from a file.

    Args:
        file_path (str): The full path to the binary file.

    Returns:
        dict: A dictionary containing the extracted metadata.
              Example: {'type': 'BINARY', 'name': 'example', 'ext': 'png', 'size': 12345}
    """
    file_name_with_ext = os.path.basename(file_path)
    file_name, file_ext = os.path.splitext(file_name_with_ext)
    file_ext = file_ext.lstrip('.') # Remove leading dot from extension

    file_size = os.path.getsize(file_path)

    metadata = {
        'type': 'BINARY',
        'name': file_name_with_ext, # Storing full name with extension as requested for FASTA header
        'ext': file_ext,
        'size': file_size
    }
    return metadata

def serialize_metadata_to_fasta_header(metadata):
    """
    Converts a metadata dictionary into a string suitable for a FASTA header.

    Args:
        metadata (dict): A dictionary containing metadata (e.g., {'type': 'BINARY', 'name': 'example.png', 'ext': 'png', 'size': 12345}).

    Returns:
        str: A string formatted as a FASTA header, e.g., '>BINARY;name=example.png;ext=png;size=12345'.
             The leading '>' is included.
    """
    header_parts = []
    # The first part is always the type, which is mandatory
    if 'type' in metadata:
        header_parts.append(metadata['type'])
    else:
        # If type is missing, this indicates malformed metadata, return a default/error header
        return '>UNKNOWN_TYPE'

    # Add other metadata fields as key=value pairs
    for key, value in metadata.items():
        if key != 'type' and value is not None:
            # Ensure values are string-convertible
            header_parts.append(f"{key}={value}")

    # Join with semicolon and prepend '>' to form the final FASTA header string
    return '>' + ';'.join(header_parts)

def parse_metadata_from_fasta_header(header_string):
    """
    Parses a FASTA header string to reconstruct the original metadata dictionary.

    Args:
        header_string (str): The FASTA header string, e.g., '>BINARY;name=example.png;ext=png;size=12345'.

    Returns:
        dict: A dictionary containing the parsed metadata. Returns an empty dict
              or a dict with an 'error' key if parsing fails or is incomplete.
    """
    metadata = {}

    # Remove leading '>' if present
    if header_string.startswith('>'):
        header_string = header_string[1:]

    parts = header_string.split(';')

    if not parts:
        metadata['error'] = 'Empty header string provided.'
        return metadata

    # The first part is assumed to be the 'type'
    if parts[0]: # Ensure it's not an empty string
        metadata['type'] = parts[0]
    else:
        metadata['error'] = 'Metadata type is missing or malformed.'

    # Parse other key=value pairs
    for part in parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            # Attempt to convert known numeric fields
            if key == 'size':
                try:
                    metadata[key] = int(value)
                except ValueError:
                    metadata[key] = value # Keep as string if not a valid int
            else:
                metadata[key] = value
        # Handle cases where a part might be empty or malformed (e.g., 'key=')
        elif part.strip(): # If it's not just whitespace, might be malformed but worth noting
            if 'malformed_parts' not in metadata:
                metadata['malformed_parts'] = []
            metadata['malformed_parts'].append(part.strip())

    return metadata

# =========================
# FASTA FUNCTIONS
# =========================

def generate_fasta_string(seq, header, line_width=60):
    """
    Generates a FASTA formatted string from a nucleotide sequence.
    Accepts a pre-formatted header string (including '>').
    """
    # Ensure header starts with '>', if not, add it
    if not header.startswith('>'):
        header = '>' + header
    out_lines = [header]
    for i in range(0, len(seq), line_width):
        out_lines.append("".join(seq[i:i+line_width]))
    return "\n".join(out_lines) + "\n"

# =========================
# DUMMY VISUALIZATION FUNCTION
# =========================

def visualize_nucleotide_sequence(seq, title="AGC-128 Sequence", checksum_length=0, error_index=-1):
    """
    Плейсхолдър – няма графика, само показва информация.
    """
    info_message = f"Title: {title}\n"
    info_message += f"Sequence Length: {len(seq)} nucleotides\n"
    if checksum_length > 0:
        info_message += f"Checksum Length: {checksum_length} nucleotides\n"
        info_message += f"Checksum Nucleotides: {' '.join(seq[-checksum_length:])}\n"
    if error_index != -1:
        info_message += f"Highlighted Error at index: {error_index} (nucleotide: {seq[error_index]})\n"
    info_message += (
        "\n(Visualization functionality is a placeholder in this environment. "
        "Run locally for full matplotlib visualization.)"
    )

    messagebox.showinfo(
        "Visualize Sequence (Placeholder)",
        info_message
    )

# =========================
# SELF-RECOVERY FUNCTIONS
# =========================

def nucleotide_sequence_to_binary_string(nucleotide_sequence):
    """
    Converts a list of nucleotides (e.g., ['A', 'T', 'G', 'C']) into a concatenated binary string.
    Example: ['G', 'C'] -> '1100'
    """
    binary_chunks = []
    for nuc in nucleotide_sequence:
        int_val = nuc_to_int.get(nuc)
        if int_val is None:
            raise ValueError(f"Invalid nucleotide character encountered: {nuc}")
        binary_chunks.append(f"{int_val:02b}")
    return "".join(binary_chunks)


def detect_no_triple_rule_violations(nucleotide_sequence):
    """
    Detects violations of the 'No-Triple Rule' (111 or 000) in a nucleotide sequence.

    Args:
        nucleotide_sequence (list): A list of nucleotide characters (e.g., ['A', 'T', 'G', 'C']).

    Returns:
        list: A list of dictionaries, where each dictionary contains the 'binary_index',
              'nucleotide_index', and 'violation_type' for each detected violation.
    """
    binary_string = nucleotide_sequence_to_binary_string(nucleotide_sequence)
    violations = []

    for i in range(len(binary_string) - 2):
        segment = binary_string[i : i + 3]
        if segment == "000":
            violations.append({
                "binary_index": i,
                "nucleotide_index": i // 2, # Each nucleotide is 2 bits
                "violation_type": "000_triple"
            })
        elif segment == "111":
            violations.append({
                "binary_index": i,
                "nucleotide_index": i // 2,
                "violation_type": "111_triple"
            })

    return violations


def detect_deterministic_next_bit_rule_violations(nucleotide_sequence):
    """
    Detects violations of the 'Deterministic-Next-Bit Rule' in a nucleotide sequence.
    - After '11' -> the next bit must be '0'
    - After '00' -> the next bit must be '1'

    Args:
        nucleotide_sequence (list): A list of nucleotide characters (e.g., ['A', 'T', 'G', 'C']).

    Returns:
        list: A list of dictionaries, where each dictionary contains the 'binary_index',
              'nucleotide_index', and 'violation_type' for each detected violation.
    """
    binary_string = nucleotide_sequence_to_binary_string(nucleotide_sequence)
    violations = []

    for i in range(len(binary_string) - 2):
        prefix = binary_string[i : i + 2]
        next_bit = binary_string[i + 2]

        if prefix == "11" and next_bit == "1":
            violations.append({
                "binary_index": i,
                "nucleotide_index": i // 2,
                "violation_type": "11_followed_by_1_expected_0"
            })
        elif prefix == "00" and next_bit == "0":
            violations.append({
                "binary_index": i,
                "nucleotide_index": i // 2,
                "violation_type": "00_followed_by_0_expected_1"
            })

    return violations


def detect_sum_2_rule_violations(nucleotide_sequence):
    """
    Detects violations of the 'Sum-2 Rule' by checking for invalid nucleotide characters.
    In this context, 'Sum-2 Rule' is interpreted as ensuring all nucleotides are valid ('A', 'T', 'G', 'C').

    Args:
        nucleotide_sequence (list): A list of nucleotide characters (e.g., ['A', 'T', 'G', 'C']).

    Returns:
        list: A list of dictionaries, where each dictionary contains the 'nucleotide_index'
              and 'violation_type' for each detected invalid character.
    """
    violations = []
    valid_nucleotides = {'A', 'T', 'G', 'C'}

    for i, nuc in enumerate(nucleotide_sequence):
        if nuc not in valid_nucleotides:
            violations.append({
                "nucleotide_index": i,
                "violation_type": "invalid_nucleotide_character"
            })

    return violations


def find_corrupted_regions(nucleotide_sequence, has_checksum=False):
    """
    Combines all error detection rules (Sum-2, No-Triple, Deterministic-Next-Bit)
    and checksum verification to identify corrupted regions in a nucleotide sequence.

    Args:
        nucleotide_sequence (list): The list of nucleotide characters.
        has_checksum (bool): True if a 2-nucleotide checksum is expected at the end
                             of the sequence, False otherwise.

    Returns:
        list: A list of dictionaries, each describing a detected violation or corrupted region.
              Each dictionary will include at least 'type' and 'indices' (list of int).
    """
    all_violations = []
    unique_violation_signatures = set() # To avoid duplicate violation entries (e.g., same type at same index)

    # 1. Check Sum-2 Rule (invalid nucleotide characters)
    sum_2_violations = detect_sum_2_rule_violations(nucleotide_sequence)
    for v in sum_2_violations:
        signature = (v['violation_type'], v['nucleotide_index'])
        if signature not in unique_violation_signatures:
            all_violations.append({'type': v['violation_type'], 'indices': [v['nucleotide_index']]})
            unique_violation_signatures.add(signature)

    # 2. Check No-Triple Rule (000 or 111 in binary string)
    no_triple_violations = detect_no_triple_rule_violations(nucleotide_sequence)
    for v in no_triple_violations:
        # No-Triple rule violations are often overlapping; represent by the start nucleotide index
        signature = (v['violation_type'], v['nucleotide_index'])
        if signature not in unique_violation_signatures:
            all_violations.append({'type': v['violation_type'], 'indices': [v['nucleotide_index']], 'binary_start': v['binary_index']})
            unique_violation_signatures.add(signature)

    # 3. Check Deterministic-Next-Bit Rule
    deterministic_next_violations = detect_deterministic_next_bit_rule_violations(nucleotide_sequence)
    for v in deterministic_next_violations:
        # Deterministic-Next-Bit violations cover a 3-bit window starting at binary_index
        signature = (v['violation_type'], v['nucleotide_index'])
        if signature not in unique_violation_signatures:
            all_violations.append({'type': v['violation_type'], 'indices': [v['nucleotide_index']], 'binary_start': v['binary_index']})
            unique_violation_signatures.add(signature)

    # 4. Check Checksum if expected
    if has_checksum:
        # Assume checksum is the last two nucleotides if present
        if len(nucleotide_sequence) < 2:
            signature = ('checksum_error', 'sequence_too_short_for_checksum')
            if signature not in unique_violation_signatures:
                all_violations.append({'type': 'checksum_error', 'description': 'Sequence too short to contain a checksum.', 'indices': []})
                unique_violation_signatures.add(signature)
        else:
            is_checksum_valid = verify_genetic_checksum(nucleotide_sequence)
            if not is_checksum_valid:
                checksum_start_index = len(nucleotide_sequence) - 2
                checksum_indices = [checksum_start_index, checksum_start_index + 1]
                signature = ('checksum_invalid', tuple(checksum_indices))
                if signature not in unique_violation_signatures:
                    all_violations.append({'type': 'checksum_invalid', 'indices': checksum_indices, 'description': 'Checksum does not match data.'})
                    unique_violation_signatures.add(signature)

    return all_violations


def binary_string_to_nucleotide_sequence(binary_string):
    """
    Converts a binary string (e.g., '1100') back into a list of nucleotides (e.g., ['G', 'C']).
    Each 2 bits are converted to one nucleotide.
    """
    if len(binary_string) % 2 != 0:
        raise ValueError("Binary string length must be a multiple of 2 to convert to nucleotides.")

    nucleotide_sequence = []
    for i in range(0, len(binary_string), 2):
        two_bit_chunk = binary_string[i:i+2]
        # Convert 2-bit string to integer, then to nucleotide
        int_val = int(two_bit_chunk, 2)
        nucleotide = int_to_nuc[int_val]
        nucleotide_sequence.append(nucleotide)
    return nucleotide_sequence


def deterministic_reconstruction_segment(nucleotide_sequence, violation):
    """
    Attempts to deterministically reconstruct a small segment of the sequence
    based on 'Deterministic-Next-Bit Rule' violations.

    Args:
        nucleotide_sequence (list): The original list of nucleotide characters.
        violation (dict): A dictionary describing the violation, typically from
                          `detect_deterministic_next_bit_rule_violations`.
                          Expected keys: 'violation_type', 'binary_index'.

    Returns:
        list: The reconstructed nucleotide sequence if a correction was applied,
              otherwise the original sequence.
        bool: True if a correction was applied, False otherwise.
    """
    binary_string_list = list(nucleotide_sequence_to_binary_string(nucleotide_sequence))
    binary_index = violation['binary_index']
    violation_type = violation['violation_type']
    corrected = False

    # Ensure binary_index + 2 is within bounds for bit flipping
    if binary_index + 2 < len(binary_string_list):
        if violation_type == "11_followed_by_1_expected_0":
            # If '11' is followed by '1', it should be '0'. Flip '1' to '0'.
            if binary_string_list[binary_index + 2] == '1':
                binary_string_list[binary_index + 2] = '0'
                corrected = True
        elif violation_type == "00_followed_by_0_expected_1":
            # If '00' is followed by '0', it should be '1'. Flip '0' to '1'.
            if binary_string_list[binary_index + 2] == '0':
                binary_string_list[binary_index + 2] = '1'
                corrected = True

    if corrected:
        reconstructed_sequence = binary_string_to_nucleotide_sequence("".join(binary_string_list))
        return reconstructed_sequence, corrected
    else:
        return nucleotide_sequence, corrected

def generate_candidate_variants(corrupted_segment, segment_start_index, full_original_sequence, max_variants=100, max_depth=None):
    """
    (Conceptual) Generates all structurally valid alternative nucleotide sequences
    for a larger corrupted segment using recursive backtracking on a bit level.
    This is a placeholder function reflecting the design discussed in the README.

    Args:
        corrupted_segment (list): The original nucleotide segment identified as corrupted.
        segment_start_index (int): The starting index of the corrupted_segment in the full_original_sequence.
        full_original_sequence (list): The entire original nucleotide sequence.
        max_variants (int): Maximum number of variants to generate.
        max_depth (int): Maximum number of nucleotides to change/generate.

    Returns:
        list: A list of lists of strings, where each inner list is a valid candidate
              nucleotide sequence for the segment.
    """
    # Placeholder: In a full implementation, this would involve a complex
    # recursive backtracking algorithm validating against No-Triple and
    # Deterministic-Next-Bit rules at each step, considering context from
    # full_original_sequence at the segment boundaries.

    # For now, return a single variant (the original corrupted segment) or an empty list
    # if no complex recovery logic is needed/possible.
    print(f"\n[Conceptual Function Call] generate_candidate_variants called for segment starting at index {segment_start_index} with length {len(corrupted_segment)}")
    print("This function is conceptual and would perform a backtracking search for valid variations.")

    # Example: If there was a simple single-nucleotide flip to consider (illustrative)
    if len(corrupted_segment) == 1 and corrupted_segment[0] == 'X': # 'X' representing an invalid nuc
        return [['A'], ['T'], ['G'], ['C']] # Try all possibilities

    return [corrupted_segment] # Return the original segment as the only 'variant' for now

def select_best_variant_with_checksum(candidate_variants, original_sequence_prefix, original_sequence_suffix):
    """
    Selects the best candidate variant based on the validity of the genetic checksum
    of the full reconstructed sequence.

    Args:
        candidate_variants (list): A list of nucleotide sequences (lists of chars),
                                   each representing a potential correction for a segment.
        original_sequence_prefix (list): The part of the original sequence before the corrupted segment.
        original_sequence_suffix (list): The part of the original sequence after the corrupted segment.
                                   This *must* include the original 2-nucleotide checksum if one was present.

    Returns:
        list or None: The nucleotide sequence of the best variant that results in a valid checksum,
                      or None if no variant yields a valid checksum.
    """
    best_variant = None

    for candidate_segment in candidate_variants:
        # Construct the full reconstructed sequence by integrating the candidate segment
        full_reconstructed_sequence = original_sequence_prefix + candidate_segment + original_sequence_suffix

        # Verify the genetic checksum of the full sequence
        # The `verify_genetic_checksum` function assumes the checksum is the last two nucleotides
        # of the *entire* sequence passed to it.
        # Therefore, `original_sequence_suffix` must include these two checksum nucleotides
        # if a checksum was originally present for the full sequence.
        if len(full_reconstructed_sequence) >= 2:
            is_checksum_valid = verify_genetic_checksum(full_reconstructed_sequence)

            if is_checksum_valid:
                # If a valid checksum is found, this is the best variant, return immediately.
                best_variant = candidate_segment
                return best_variant

    # If no variant resulted in a valid checksum, return None
    return None

def attempt_self_recovery(nucleotide_sequence, has_checksum=False):
    """
    Orchestrates the self-recovery process for a nucleotide sequence.

    Args:
        nucleotide_sequence (list): The original list of nucleotide characters to recover.
        has_checksum (bool): True if a 2-nucleotide checksum is expected at the end
                             of the sequence, False otherwise.

    Returns:
        tuple: A tuple containing:
            - list: The recovered nucleotide sequence.
            - dict: A recovery report including fixed errors, unresolved issues, and checksum status.
    """
    recovered_sequence = copy.deepcopy(nucleotide_sequence)
    fixed_errors = 0
    unresolved_issues = []

    # Find initial corrupted regions
    initial_violations = find_corrupted_regions(recovered_sequence, has_checksum=has_checksum)

    # Iterate through each violation and attempt deterministic correction
    for violation in initial_violations:
        # Deterministic-Next-Bit Rule violations can be fixed deterministically
        if violation['type'] in ["11_followed_by_1_expected_0", "00_followed_by_0_expected_1"]:
            # Attempt deterministic correction
            temp_sequence, corrected = deterministic_reconstruction_segment(recovered_sequence, violation)
            if corrected:
                recovered_sequence = temp_sequence
                fixed_errors += 1
            else:
                # If deterministic reconstruction couldn't apply a fix (e.g., already fixed, or not simple enough)
                unresolved_issues.append({
                    'original_violation': violation,
                    'attempted_fix': 'deterministic_failed',
                    'reason': 'Deterministic reconstruction did not apply a fix.'
                })
        elif violation['type'] == 'invalid_nucleotide_character':
            # This type of error is not deterministically fixable without more context
            # or candidate generation logic. For now, mark as unresolved.
            unresolved_issues.append({
                'original_violation': violation,
                'attempted_fix': 'none',
                'reason': 'Invalid nucleotide character, requires more complex recovery logic (e.g., generate_candidate_variants).'
            })
        elif violation['type'] in ['000_triple', '111_triple']:
            # Triple violations are also not deterministically fixable with a single bit flip
            unresolved_issues.append({
                'original_violation': violation,
                'attempted_fix': 'none',
                'reason': 'No-Triple Rule violation, requires more complex recovery logic (e.g., generate_candidate_variants).'
            })
        elif violation['type'] in ['checksum_invalid', 'checksum_error']:
            # Checksum issues are handled separately in the final report, not as individual "unresolved issues" during iterative fixing
            pass # These will be captured by final_checksum_status
        else:
            # Catch-all for any other unforeseen violation types
            unresolved_issues.append({
                'original_violation': violation,
                'attempted_fix': 'none',
                'reason': 'Unhandled violation type.'
            })

    # After attempting deterministic fixes, re-evaluate all violations to get the final state
    final_violations_after_deterministic = find_corrupted_regions(recovered_sequence, has_checksum=has_checksum)

    # Update unresolved_issues based on persistent or newly revealed problems after deterministic pass
    # This part is simplified: if an issue (identified by its type and relevant indices) is still present,
    # and it hasn't already been explicitly marked as 'deterministic_failed', add it to unresolved.
    # A more sophisticated approach would track changes to each specific violation.
    for f_violation in final_violations_after_deterministic:
        # Check if this violation was already handled (e.g. fixed or marked as deterministic_failed)
        is_resolved_or_tracked = False
        if f_violation['type'] in ["11_followed_by_1_expected_0", "00_followed_by_0_expected_1"]:
            # If it's a deterministic type, and we tried to fix it, it should have been corrected
            # or already marked 'deterministic_failed'. If it's still here, it's a persistent problem.
            # This simple check assumes deterministic fixes don't create *new* identical violations
            pass # Will be implicitly part of unresolved_issues if not fixed by deterministic_reconstruction_segment
        else:
            # For other types (Sum-2, No-Triple, etc.), if still present, it's unresolved by simple means
            # Avoid adding duplicates if already marked during initial pass
            if not any(u_issue.get('original_violation') == f_violation for u_issue in unresolved_issues):
                unresolved_issues.append({
                    'original_violation': f_violation,
                    'attempted_fix': 'none_after_deterministic_pass',
                    'reason': 'Violation persisted after deterministic attempts or was not fixable deterministically.'
                })

    # Construct recovery report
    final_checksum_status = None
    if has_checksum:
        # The verify_genetic_checksum function will automatically extract the data part
        final_checksum_status = verify_genetic_checksum(recovered_sequence)

    recovery_report = {
        'fixed_errors_count': fixed_errors,
        'unresolved_issues': unresolved_issues,
        'final_checksum_status': final_checksum_status,
        'recovered_sequence_length': len(recovered_sequence)
    }

    return recovered_sequence, recovery_report

def calculate_recovery_confidence(recovery_report):
    """
    Calculates a confidence score for the self-recovery process (0-100).

    Args:
        recovery_report (dict): The report generated by attempt_self_recovery.
                                Expected keys: 'fixed_errors_count', 'unresolved_issues',
                                'final_checksum_status'.

    Returns:
        int: The recovery confidence score, an integer between 0 and 100.
    """
    confidence_score = 0

    # 1. Base score based on final checksum status (most critical factor)
    if recovery_report['final_checksum_status'] is True:
        confidence_score = 100  # Highest confidence if checksum is valid
    elif recovery_report['final_checksum_status'] is False:
        confidence_score = 0   # Lowest confidence if checksum is invalid
    else: # final_checksum_status is None (checksum not used/expected)
        confidence_score = 50  # Neutral starting point

    # 2. Adjust for fixed errors (only if checksum is not explicitly false, to avoid negative scores from start)
    if recovery_report['final_checksum_status'] is not False:
        # Each fixed error adds a small positive contribution, up to a limit
        confidence_score += min(recovery_report['fixed_errors_count'] * 2, 10) # Max +10 points for fixed errors

    # 3. Penalize for unresolved issues
    # Penalties are more severe if checksum is not valid or not used.
    for issue in recovery_report['unresolved_issues']:
        violation_type = issue.get('original_violation', {}).get('type')
        # attempted_fix_status = issue.get('attempted_fix') # Not directly used for scoring here, but useful for debugging

        if recovery_report['final_checksum_status'] is True:
            # Minor penalties if checksum is true, as it implies overall integrity was restored
            if violation_type == 'invalid_nucleotide_character':
                confidence_score -= 3
            elif violation_type in ['000_triple', '111_triple']:
                confidence_score -= 2
            elif violation_type in ['11_followed_by_1_expected_0', '00_followed_by_0_expected_1']:
                confidence_score -= 1
            # For any other lingering issues (that checksum still validates over), very minor penalty
            else:
                confidence_score -= 1
        elif recovery_report['final_checksum_status'] is None:
            # More significant penalties if no checksum validation to fall back on
            if violation_type == 'invalid_nucleotide_character':
                confidence_score -= 15
            elif violation_type in ['000_triple', '111_triple']:  
                confidence_score -= 10
            elif violation_type in ['11_followed_by_1_expected_0', '00_followed_by_0_expected_1']:
                confidence_score -= 7
            # Checksum issues themselves should not appear here if final_checksum_status is None
            # but if they do, apply heavy penalty
            elif 'checksum' in violation_type:
                confidence_score -= 20
            else:
                confidence_score -= 5
        # If final_checksum_status is False, score is already 0, no further penalties needed.

    # 4. Clamp the score between 0 and 100
    confidence_score = max(0, min(100, confidence_score))

    return confidence_score


# =========================
# AGC-256 CONCEPTUAL FUNCTIONS/CLASSES
# =========================

# Law I: Windowing Functions
def calculate_2_bit_windows(bit_stream_length):
    """
    Calculates the number of overlapping 2-bit windows for a given bit stream length.
    Formula: bit_stream_length - 2 + 1
    """
    if bit_stream_length < 2:
        return 0
    return bit_stream_length - 2 + 1

def calculate_3_bit_windows(bit_stream_length):
    """
    Calculates the number of overlapping 3-bit windows for a given bit stream length.
    Formula: bit_stream_length - 3 + 1
    """
    if bit_stream_length < 3:
        return 0
    return bit_stream_length - 3 + 1

def calculate_4_bit_windows(bit_stream_length):
    """
    Calculates the number of overlapping 4-bit windows for a given bit stream length.
    Formula: bit_stream_length - 4 + 1
    """
    if bit_stream_length < 4:
        return 0
    return bit_stream_length - 4 + 1

# Law II: Nucleotide Connection Functions
def calculate_ordered_nucleotide_pairs(num_nucleotides):
    """
    Calculates the number of ordered nucleotide pairs for a given number of nucleotides.
    Formula: (n * (n - 1)) / 2
    """
    if num_nucleotides < 2:
        return 0
    return (num_nucleotides * (num_nucleotides - 1)) / 2

def calculate_connections_with_2_modes(num_nucleotides):
    """
    Calculates the number of connections with 2 modes for a given number of nucleotides.
    Formula: num_nucleotides * (num_nucleotides - 1)
    """
    if num_nucleotides < 1:
        return 0
    return num_nucleotides * (num_nucleotides - 1)

# Law III: Fractal Container Data Structure
class FractalCube:
    """
    Conceptual Python class to represent the nested structure of an AGC-256 'Cube',
    modeling a core, internal context, and external context.
    """
    def __init__(self, core=None, internal_context=None, external_context=None):
        """
        Initializes the FractalCube with its conceptual components.

        Args:
            core: Represents the central 2-bit motif (e.g., 'AG', 'CT').
            internal_context: Represents the two 1-bit motifs directly surrounding the core (e.g., ('0', '1')).
            external_context: Represents the two 1-bit motifs further out,
                              providing external context to the internal_context (e.g., ('0', '1')).
        """
        self.core = core if core is not None else ""
        self.internal_context = internal_context if internal_context is not None else ("", "")
        self.external_context = external_context if external_context is not None else ("", "")

    def __repr__(self):
        """
        Provides a human-readable representation of the FractalCube structure.
        """
        return (
            f"FractalCube("
            f"  core='{self.core}',\n"
            f"  internal_context=('{self.internal_context[0]}', '{self.internal_context[1]}'),\n"
            f"  external_context=('{self.external_context[0]}', '{self.external_context[1]}')"
            f")"
        )

# Law IV: Fractal Function
def calculate_fractal_motives(bit_stream_length, motive_size):
    """
    Calculates the number of motives of size 'motive_size' (k)
    from a bit stream of length 'bit_stream_length' (n).
    Formula: (n - k + 1) * (2 ** k)
    """
    if bit_stream_length < motive_size:
        return 0
    return (bit_stream_length - motive_size + 1) * (2 ** motive_size)

# Law V: AGC-256 Cube Derivation
def simulate_agc256_cube_derivation(core_nucleotide, left_context_bit, right_context_bit):
    """
    Simulates how a 4-bit AGC-256 cube is derived from a 2-bit core and 1-bit context bits.

    Args:
        core_nucleotide (str): A 2-bit nucleotide ('A', 'T', 'G', 'C').
        left_context_bit (str): A 1-bit ('0' or '1'). This will be the first bit of the 4-bit cube.
        right_context_bit (str): A 1-bit ('0' or '1'). This will be the last bit of the 4-bit cube.

    Returns:
        FractalCube: A FractalCube object representing the derived 4-bit structure.
                     For a 4-bit cube, internal and external contexts will be identical
                     as they represent the flanking bits of the core.
    """
    # 1. Convert core_nucleotide to its 2-bit binary string
    core_int_val = nuc_to_int.get(core_nucleotide)
    if core_int_val is None:
        raise ValueError(f"Invalid core nucleotide: {core_nucleotide}")
    core_binary = f"{core_int_val:02b}"

    # 2. Combine bits into a single 4-bit binary string
    # The full 4-bit binary string is (left_context_bit) + (core_binary_0) + (core_binary_1) + (right_context_bit)
    four_bit_binary_string = left_context_bit + core_binary[0] + core_binary[1] + right_context_bit

    # 3. Create a FractalCube object
    # For a 4-bit cube, the 'internal_context' (flanking the core) and 'external_context'
    # (outermost bits) are the same, as there's no layer 'further out'.
    fractal_cube = FractalCube(
        core=core_nucleotide,
        internal_context=(left_context_bit, right_context_bit),
        external_context=(left_context_bit, right_context_bit) # Explicitly set to be the same as internal for 4-bit cube
    )

    return fractal_cube


# =========================
# GUI SETUP
# =========================

def setup_gui():
    global current_encoded_nucleotide_sequence

    root = tk.Tk()
    root.title("AGC-128 Notepad")
    root.geometry("1050x600") # Set initial window size

    # Frame for encoding version selection
    version_frame = tk.Frame(root)
    version_frame.pack(pady=5, anchor='w')

    tk.Label(version_frame, text="Encoding/Decoding Version:").pack(side=tk.LEFT)
    version_var = tk.StringVar(value="v1_ascii")  # Default to v1 (ASCII)

    v1_radio = tk.Radiobutton(version_frame, text="v1 (ASCII)", variable=version_var, value="v1_ascii")
    v1_radio.pack(side=tk.LEFT, padx=5)

    v2_radio = tk.Radiobutton(version_frame, text="v2 (Unicode)", variable=version_var, value="v2_unicode")
    v2_radio.pack(side=tk.LEFT, padx=5)

    # Configure text_widget with undo/redo history and scrollbar
    text_frame = tk.Frame(root) # New frame for text widget and scrollbar
    text_frame.pack(expand=True, fill='both')

    text_widget = tk.Text(text_frame, wrap='word', undo=True, autoseparators=True)
    text_widget.pack(side=tk.LEFT, expand=True, fill='both')

    # Add a scrollbar
    scrollbar = tk.Scrollbar(text_frame, command=text_widget.yview)
    scrollbar.pack(side=tk.RIGHT, fill='y')
    text_widget.config(yscrollcommand=scrollbar.set)

    menubar = tk.Menu(root)
    root.config(menu=menubar)

    # ---------- FILE ----------
    file_menu = tk.Menu(menubar, tearoff=0)
    menubar.add_cascade(label="File", menu=file_menu)

    def new_file():
        text_widget.delete("1.0", tk.END)
        current_encoded_nucleotide_sequence.clear()
        messagebox.showinfo("New File", "New file created. Editor cleared.")

    def open_file():
        global current_encoded_nucleotide_sequence
        file_path = filedialog.askopenfilename(
            filetypes=[("Text files", "*.txt"), ("All files", "*.* затем")]
        )
        if file_path:
            with open(file_path, 'r', encoding='utf-8') as file:
                content = file.read()
            text_widget.delete("1.0", tk.END)
            text_widget.insert(tk.END, content)
            current_encoded_nucleotide_sequence.clear()

    def save_file():
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.* затем")]
        )
        if file_path:
            content = text_widget.get("1.0", tk.END)
            with open(file_path, 'w', encoding='utf-8') as file:
                file.write(content)

    def save_file_as():
        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.* затем")]
        )
        if file_path:
            content = text_widget.get("1.0", tk.END)
            with open(file_path, 'w', encoding='utf-8') as file:
                file.write(content)

    file_menu.add_command(label="New", command=new_file)
    file_menu.add_command(label="Open", command=open_file)
    file_menu.add_command(label="Save", command=save_file)
    file_menu.add_command(label="Save As...", command=save_file_as)
    file_menu.add_separator()
    file_menu.add_command(label="Exit", command=root.quit)

    # ---------- EDIT MENU ----------
    edit_menu = tk.Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Edit", menu=edit_menu)

    def undo_action():
        try:
            text_widget.edit_undo()
        except tk.TclError:
            pass # Cannot undo

    def redo_action():
        try:
            text_widget.edit_redo()
        except tk.TclError:
            pass # Cannot redo

    def cut_action():
        text_widget.event_generate('<<Cut>>')

    def copy_action():
        text_widget.event_generate('<<Copy>>')

    def paste_action():
        text_widget.event_generate('<<Paste>>')

    def delete_action():
        try:
            text_widget.delete(tk.SEL_FIRST, tk.SEL_LAST)
        except tk.TclError: # No text selected
            pass

    def select_all_action():
        text_widget.tag_add(tk.SEL, '1.0', tk.END)
        text_widget.mark_set(tk.INSERT, '1.0')
        text_widget.see(tk.INSERT) # Scroll to the beginning

    edit_menu.add_command(label="Undo", command=undo_action)
    edit_menu.add_command(label="Redo", command=redo_action)
    edit_menu.add_separator()
    edit_menu.add_command(label="Cut", command=cut_action)
    edit_menu.add_command(label="Copy", command=copy_action)
    edit_menu.add_command(label="Paste", command=paste_action)
    edit_menu.add_command(label="Delete", command=delete_action)
    edit_menu.add_separator()
    edit_menu.add_command(label="Select All", command=select_all_action)

    # ---------- CONTEXT MENU ----------
    def show_context_menu(event):
        context_menu = tk.Menu(text_widget, tearoff=0)
        context_menu.add_command(label="Cut", command=cut_action)
        context_menu.add_command(label="Copy", command=copy_action)
        context_menu.add_command(label="Paste", command=paste_action)
        context_menu.add_separator()
        context_menu.add_command(label="Select All", command=select_all_action)
        context_menu.add_command(label="Clear", command=lambda: text_widget.delete('1.0', tk.END))
        try:
            context_menu.tk_popup(event.x_root, event.y_root)
        finally:
            context_menu.grab_release()

    text_widget.bind("<Button-3>", show_context_menu)

    # ---------- ENCODE ----------
    encode_menu = tk.Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Encode", menu=encode_menu)

    def encode_to_fasta_action():
        global current_encoded_nucleotide_sequence

        input_text = text_widget.get("1.0", tk.END).strip()
        if not input_text:
            messagebox.showwarning("No Input", "Please enter text to encode in the editor.")
            return

        fasta_id = simpledialog.askstring("FASTA Identifier", "Enter FASTA header ID:")
        if not fasta_id:
            messagebox.showwarning("Missing ID", "FASTA identifier cannot be empty.")
            return

        add_checksum = messagebox.askyesno("Checksum Option", "Do you want to add a genetic checksum?")

        try:
            selected_version = version_var.get()
            if selected_version == "v1_ascii":
                nucleotide_sequence_temp = string_to_nucleotide_sequence_v1(input_text)
            else:  # v2_unicode
                nucleotide_sequence_temp = encode_string_to_unicode_tagc_sequence(input_text)

            if add_checksum:
                processed_sequence = add_genetic_checksum(nucleotide_sequence_temp)
            else:
                processed_sequence = nucleotide_sequence_temp

            current_encoded_nucleotide_sequence[:] = processed_sequence

            # Use the provided fasta_id as the header directly for text encoding
            fasta_output = generate_fasta_string(
                processed_sequence,
                fasta_id, # This is the raw ID, generate_fasta_string will prepend '>'
                line_width=60
            )

            save_path = filedialog.asksaveasfilename(
                defaultextension=".fasta",
                filetypes=[("FASTA files", "*.fasta"), ("All files", "*.* затем")],
                title="Save Encoded FASTA As"
            )
            if save_path:
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(fasta_output)
                messagebox.showinfo("Success", f"FASTA encoded and saved to {save_path}")
            else:
                messagebox.showinfo("Cancelled", "FASTA save operation cancelled.")
        except Exception as e:
            messagebox.showerror("Encoding Error", f"An error occurred during encoding: {e}")

    encode_menu.add_command(label="Encode Text to AGC-128 FASTA", command=encode_to_fasta_action) # Renamed for clarity

    # Binary Submenu under Encode
    binary_encode_menu = tk.Menu(encode_menu, tearoff=0)
    encode_menu.add_cascade(label="Binary", menu=binary_encode_menu)

    def encode_binary_file_action():
        global current_encoded_nucleotide_sequence

        file_path = filedialog.askopenfilename(
            title="Select Binary File to Encode",
            filetypes=[("All files", "*.* затем")]
        )
        if not file_path:
            messagebox.showinfo("Cancelled", "Binary file encoding cancelled.")
            return

        add_checksum = messagebox.askyesno("Checksum Option", "Do you want to add a genetic checksum?")

        try:
            # a. Read file content as raw bytes
            with open(file_path, 'rb') as f:
                raw_bytes_content = f.read()

            # b. Call extract_file_metadata()
            metadata = extract_file_metadata(file_path)

            # c. Call bytes_to_nucleotide_sequence() to convert raw binary content
            nucleotide_sequence_temp = bytes_to_nucleotide_sequence(raw_bytes_content)

            if add_checksum:
                processed_sequence = add_genetic_checksum(nucleotide_sequence_temp)
            else:
                processed_sequence = nucleotide_sequence_temp

            current_encoded_nucleotide_sequence[:] = processed_sequence

            # d. Call serialize_metadata_to_fasta_header() to create FASTA header
            fasta_header = serialize_metadata_to_fasta_header(metadata)

            # e. Use generate_fasta_string() to construct final FASTA content
            fasta_output = generate_fasta_string(
                processed_sequence,
                fasta_header, # Use the pre-formatted metadata header
                line_width=60
            )

            # f. Prompt user for save location and write
            save_path = filedialog.asksaveasfilename(
                defaultextension=".fasta",
                filetypes=[("FASTA files", "*.fasta"), ("All files", "*.* затем")],
                title="Save Encoded Binary FASTA As"
            )
            if save_path:
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(fasta_output)
                messagebox.showinfo("Success", f"Binary file encoded to FASTA and saved to {save_path}")
            else:
                messagebox.showinfo("Cancelled", "FASTA save operation cancelled.")

        except Exception as e:
            messagebox.showerror("Encoding Error", f"An error occurred during binary encoding: {e}")

    binary_encode_menu.add_command(label="Encode Binary File to AGC-128 FASTA", command=encode_binary_file_action)

    # ---------- DECODE ----------
    decode_menu = tk.Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Decode", menu=decode_menu)

    def load_and_decode_fasta_action(expected_type=None):
        global current_encoded_nucleotide_sequence

        file_path = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.* затем")]
        )
        if not file_path:
            messagebox.showinfo("Cancelled", "FASTA load operation cancelled.")
            return

        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                content = file.read()

            lines = content.splitlines()
            if not lines:
                messagebox.showwarning(
                    "Invalid FASTA",
                    "Selected file is empty or does not appear to be a valid FASTA format."
                )
                return

            fasta_header_line = lines[0]
            if not fasta_header_line.startswith('>'):
                messagebox.showwarning(
                    "Invalid FASTA",
                    "Selected file does not appear to be a valid FASTA format (missing header)."
                )
                return

            # Extract metadata from header line
            metadata = parse_metadata_from_fasta_header(fasta_header_line)
            file_type = metadata.get('type', 'TEXT') # Default to TEXT if type is not specified
            original_filename = metadata.get('name', 'decoded_output.txt') # Default filename for saving binary

            # Check if expected_type matches actual file_type from metadata
            if expected_type and file_type != expected_type:
                # Display a warning but proceed based on actual metadata type
                messagebox.showwarning(
                    "Type Mismatch",
                    f"Expected a {expected_type} FASTA file, but found type '{file_type}' in metadata.\n"
                    "Attempting to decode as {file_type} anyway."
                )

            # Extract sequence, ignore header(s), keep only A/T/G/C
            seq_raw = "".join(line.strip() for line in lines[1:] if not line.startswith(">"))
            valid = {'A', 'T', 'G', 'C'}
            extracted_nucs_list = [c for c in seq_raw if c in valid]

            if not extracted_nucs_list:
                messagebox.showwarning("Empty Sequence", "No nucleotide sequence found in the FASTA file.")
                return

            current_encoded_nucleotide_sequence[:] = extracted_nucs_list

            sequence_to_decode = list(extracted_nucs_list) # Use a copy to allow modification
            checksum_info = ""

            # --- MODIFIED CHECKSUM HANDLING ---
            ask_if_checksum_present = messagebox.askyesno(
                "Checksum Query",
                "Is a 2-nucleotide genetic checksum expected at the end of this sequence?"
            )

            if ask_if_checksum_present:
                if len(extracted_nucs_list) < 2:
                    messagebox.showwarning("Checksum Error", "Sequence is too short to contain a 2-nucleotide checksum.")
                else:
                    is_valid_checksum = verify_genetic_checksum(extracted_nucs_list)
                    checksum_info = f"\nChecksum valid: {is_valid_checksum}"
                    if is_valid_checksum:
                        messagebox.showinfo("Checksum Status", f"Checksum is valid!{checksum_info}")
                        sequence_to_decode = extracted_nucs_list[:-2] # Remove checksum for decoding
                    else:
                        # If checksum is invalid, still remove it, especially for binary files,
                        # to allow further processing but warn the user.
                        messagebox.showwarning(
                            "Checksum Status",
                            f"Checksum is INVALID! Data may be corrupted.{checksum_info}\n"
                            "The checksum will be removed for decoding to allow file reconstruction, but data integrity is compromised."
                        )
                        sequence_to_decode = extracted_nucs_list[:-2] # Always remove if expected, even if invalid
            # --- END MODIFIED CHECKSUM HANDLING ---

            # Determine decoding method based on actual file_type from metadata
            if file_type == 'BINARY':
                try:
                    if len(sequence_to_decode) % 4 != 0:
                        # This warning should ideally not happen after checksum removal if original was valid
                        # but might if initial sequence was malformed.
                        messagebox.showwarning(
                            "Sequence Length Mismatch (Binary)",
                            "The binary nucleotide sequence length is not a multiple of 4.\n"
                            "Decoding might result in an incomplete last byte."
                        )

                    decoded_bytes = nucleotide_sequence_to_bytes(sequence_to_decode)
                    save_binary_path = filedialog.asksaveasfilename(
                        defaultextension=f".{metadata.get('ext', '')}",
                        initialfile=original_filename,
                        title="Save Decoded Binary File As"
                    )
                    if save_binary_path:
                        with open(save_binary_path, 'wb') as f:
                            f.write(decoded_bytes)
                        messagebox.showinfo("Decoding Success", f"Binary file successfully decoded and saved to {save_binary_path}!{checksum_info}")
                    else:
                        messagebox.showinfo("Cancelled", "Binary file save operation cancelled.")
                except Exception as e:
                    messagebox.showerror("Binary Decoding Error", f"An error occurred during binary decoding: {e}")
            else: # Assume TEXT, either v1 or v2 (TEXT is the default type for metadata handling)
                # Determine the selected version for text decoding
                selected_version = version_var.get()

                # Perform pre-decoding length check if no checksum was removed and it's V1.
                if not ask_if_checksum_present and selected_version == "v1_ascii" and len(sequence_to_decode) % 4 != 0:
                     messagebox.showwarning(
                        "Sequence Length Mismatch (V1)",
                        "The V1 ASCII nucleotide sequence length is not a multiple of 4.\n"
                        "Decoding might result in an incomplete last character."
                    )

                if selected_version == "v1_ascii":
                    decoded_text = decode_nucleotide_sequence_to_string_v1(sequence_to_decode)
                else: # v2_unicode
                    decoded_text = decode_unicode_tagc_sequence_to_string(sequence_to_decode)

                text_widget.delete("1.0", tk.END)
                text_widget.insert(tk.END, decoded_text)
                messagebox.showinfo("Decoding Success", f"FASTA file successfully loaded and decoded!{checksum_info}")

        except ValueError as ve: # Catch specific ValueError from decoding functions
            messagebox.showerror("Decoding Error (Data Integrity)", f"A data integrity error occurred during decoding: {ve}\nThis might indicate a corrupted sequence or incorrect encoding version/checksum assumption.")
        except Exception as e:
            messagebox.showerror("Decoding Error", f"An unexpected error occurred during FASTA loading or decoding: {e}")

    decode_menu.add_command(label="Load and Decode Text FASTA", command=lambda: load_and_decode_fasta_action(expected_type='TEXT')) # Renamed and added expected_type

    # Binary Submenu under Decode
    binary_decode_menu = tk.Menu(decode_menu, tearoff=0)
    decode_menu.add_cascade(label="Binary", menu=binary_decode_menu)

    binary_decode_menu.add_command(label="Decode AGC-128 to Binary File", command=lambda: load_and_decode_fasta_action(expected_type='BINARY')) # New menu item for binary decoding

    # ---------- TOOLS ----------
    tools_menu = tk.Menu(menubar, tearoff=0)
    menubar.add_cascade(label="Tools", menu=tools_menu)

    def verify_checksum_action():
        global current_encoded_nucleotide_sequence
        if not current_encoded_nucleotide_sequence:
            messagebox.showwarning("No Sequence", "No encoded nucleotide sequence is currently loaded or generated.")
            return

        # --- CHECKSUM HANDLING IN VERIFY ACTION ---
        ask_if_checksum_present = messagebox.askyesno(
            "Checksum Query",
            "Is a 2-nucleotide genetic checksum expected at the end of the current sequence?"
        )

        if ask_if_checksum_present:
            if len(current_encoded_nucleotide_sequence) < 2:
                messagebox.showwarning("Checksum Error", "The current sequence is too short to contain a 2-nucleotide checksum.")
                return

            is_valid = verify_genetic_checksum(current_encoded_nucleotide_sequence)
            messagebox.showinfo("Checksum Verification", f"Checksum valid: {is_valid}")
        else:
            messagebox.showinfo("Checksum Information", "No checksum verification performed as none was expected.")
        # --- END CHECKSUM HANDLING ---

    def visualize_action():
        global current_encoded_nucleotide_sequence
        if not current_encoded_nucleotide_sequence:
            messagebox.showwarning(
                "No Sequence",
                "No encoded nucleotide sequence is currently loaded or generated to visualize."
            )
            return

        checksum_len = 0
        sequence_for_viz = list(current_encoded_nucleotide_sequence) # Make a copy

        # --- CHECKSUM HANDLING IN VISUALIZE ACTION ---
        ask_if_checksum_present = messagebox.askyesno(
            "Checksum Query",
            "Is a 2-nucleotide genetic checksum expected at the end of the current sequence for visualization?"
        )

        if ask_if_checksum_present:
            if len(current_encoded_nucleotide_sequence) < 2:
                messagebox.showwarning("Checksum Error", "Sequence is too short to contain a 2-nucleotide checksum for visualization.")
            else:
                is_valid_checksum = verify_genetic_checksum(current_encoded_nucleotide_sequence)
                if is_valid_checksum:
                    checksum_len = 2 # Indicate to visualization to highlight last 2 nucs
                    messagebox.showinfo("Checksum Status", "Checksum is valid and will be highlighted.")
                else:
                    messagebox.showwarning("Checksum Status", "Checksum is INVALID. Will still highlight, but data may be corrupted.")
                    checksum_len = 2 # Still highlight, even if invalid
        # --- END CHECKSUM HANDLING ---

        try:
            visualize_nucleotide_sequence(
                sequence_for_viz, # Pass the original sequence, checksum_len will handle highlighting
                "Current AGC-128 Sequence",
                checksum_length=checksum_len
            )
        except Exception as e:
            messagebox.showerror("Visualization Error", f"An error occurred during visualization: {e}")

    def view_binary_representation_action():
        global current_encoded_nucleotide_sequence
        if not current_encoded_nucleotide_sequence:
            messagebox.showwarning("No Sequence", "No encoded nucleotide sequence is currently loaded or generated to view in binary.")
            return

        try:
            binary_string = nucleotide_sequence_to_binary_string(current_encoded_nucleotide_sequence)

            # Create a new Toplevel window to display the binary string
            binary_window = tk.Toplevel(root)
            binary_window.title("Binary Representation")
            binary_window.geometry("600x400")

            scrolled_text = scrolledtext.ScrolledText(binary_window, wrap='word', width=70, height=20)
            scrolled_text.pack(expand=True, fill='both', padx=10, pady=10)
            scrolled_text.insert(tk.END, binary_string)
            scrolled_text.config(state='disabled') # Make it read-only

            # Add a close button
            close_button = tk.Button(binary_window, text="Close", command=binary_window.destroy)
            close_button.pack(pady=5)

        except Exception as e:
            messagebox.showerror("Binary View Error", f"An error occurred while converting to binary: {e}")

    def attempt_recovery_action():
        global current_encoded_nucleotide_sequence
        if not current_encoded_nucleotide_sequence:
            messagebox.showwarning("No Sequence", "No encoded nucleotide sequence is currently loaded or generated for recovery.")
            return

        # Ask if checksum is present for recovery logic
        has_checksum_for_recovery = messagebox.askyesno(
            "Recovery Option",
            "Does the sequence being recovered include a 2-nucleotide checksum at the end?"
        )

        try:
            recovered_seq, recovery_report = attempt_self_recovery(current_encoded_nucleotide_sequence, has_checksum=has_checksum_for_recovery)
            confidence = calculate_recovery_confidence(recovery_report)

            # Update the global sequence if recovery was successful
            if recovered_seq != current_encoded_nucleotide_sequence:
                current_encoded_nucleotide_sequence[:] = recovered_seq
                messagebox.showinfo("Recovery Status", "Sequence was modified during recovery. The 'current_encoded_nucleotide_sequence' has been updated.")

            report_str = f"Recovery Report:\n"
            report_str += f"  Fixed Errors: {recovery_report['fixed_errors_count']}\n"
            report_str += f"  Unresolved Issues: {len(recovery_report['unresolved_issues'])} (see details below)\n"
            report_str += f"  Final Checksum Status: {recovery_report['final_checksum_status']}\n"
            report_str += f"  Recovered Sequence Length: {recovery_report['recovered_sequence_length']}\n"
            report_str += f"  Recovery Confidence: {confidence}%\n\n"
            report_str += "--- Unresolved Issues Details ---\n"
            if recovery_report['unresolved_issues']:
                for issue in recovery_report['unresolved_issues']:
                    report_str += f"    Type: {issue.get('original_violation',{}).get('type')}, "\
                                  f"Indices: {issue.get('original_violation',{}).get('indices')}, "\
                                  f"Reason: {issue.get('reason')}\n"
            else:
                report_str += "    None.\n"

            messagebox.showinfo("Self-Recovery Complete", report_str)

        except Exception as e:
            messagebox.showerror("Self-Recovery Error", f"An error occurred during self-recovery: {e}")

    tools_menu.add_command(label="Verify Checksum", command=verify_checksum_action)
    tools_menu.add_command(label="Visualize Sequence", command=visualize_action)
    tools_menu.add_command(label="View Binary Representation", command=view_binary_representation_action)
    tools_menu.add_command(label="Attempt Self-Recovery", command=attempt_recovery_action) # New recovery menu item

    root.mainloop()

# =========================
# MAIN EXECUTION BLOCK
# =========================

if __name__ == "__main__":
    # Check if running in Google Colab (or similar non-GUI environment)
    if 'google.colab' in sys.modules:
        print("Running in Google Colab environment. Tkinter GUI cannot be displayed.\n")
        print("Here's an example of how to use the core encoding/decoding functions directly:\n")

        sample_text = "Здравейте, свят!😊 123" # Updated emoji to 😊 (4-byte UTF-8)
        print(f"Original Text (V2 Unicode): {sample_text}")

        # V2 Unicode Encoding Example
        try:
            encoded_v2 = encode_string_to_unicode_tagc_sequence(sample_text)
            print(f"Encoded (V2 Unicode): {''.join(encoded_v2[:60])}{'...' if len(encoded_v2) > 60 else ''} (Total: {len(encoded_v2)} nucleotides)")

            # Add and verify checksum
            encoded_v2_with_checksum = add_genetic_checksum(encoded_v2)
            print(f"Encoded with Checksum (V2 Unicode): {''.join(encoded_v2_with_checksum[:60])}{'...' if len(encoded_v2_with_checksum) > 60 else ''} (Total: {len(encoded_v2_with_checksum)} nucleotides)")
            print(f"Checksum for V2 is valid: {verify_genetic_checksum(encoded_v2_with_checksum)}")

            decoded_v2 = decode_unicode_tagc_sequence_to_string(encoded_v2)
            print(f"Decoded (V2 Unicode): {decoded_v2}")
            print(f"V2 Encoding/Decoding successful: {sample_text == decoded_v2}")

            # AGC-256 Concepts Demonstration in Colab
            print("\n--- Demonstrating Law I: Windowing Functions ---")
            bit_len_colab = 10
            print(f"For a bit stream length of {bit_len_colab}:")
            print(f"  Number of 2-bit windows: {calculate_2_bit_windows(bit_len_colab)}")
            print(f"  Number of 3-bit windows: {calculate_3_bit_windows(bit_len_colab)}")
            print(f"  Number of 4-bit windows: {calculate_4_bit_windows(bit_len_colab)}")

            print("\n--- Demonstrating Law II: Nucleotide Connection Functions ---")
            num_nucs_colab = 4
            print(f"For {num_nucs_colab} nucleotides:")
            print(f"  Number of ordered nucleotide pairs: {calculate_ordered_nucleotide_pairs(num_nucs_colab)}")
            print(f"  Number of connections with 2 modes: {calculate_connections_with_2_modes(num_nucs_colab)}")

            print("\n--- Demonstrating Law III: Fractal Container Data Structure ---")
            print("Creating an instance of FractalCube and printing its representation:")
            # Corrected to use 1-bit strings for external_context to match definition: external_context=(four_bit_binary_string[0], four_bit_binary_string[3])
            # As per the logic in simulate_agc256_cube_derivation, for a 4-bit cube, internal and external contexts are derived from the flanking bits.
            # For a simple demo, '0' and '0' for external context bits are used.
            example_cube_colab = FractalCube(core='AG', internal_context=('0', '1'), external_context=('0', '0'))
            print(example_cube_colab)

            print("\n--- Demonstrating Law IV: Fractal Function ---")
            bit_len_fractal_colab = 10
            print(f"For a bit stream length of {bit_len_fractal_colab}:")
            print(f"  Number of motives of size 2: {calculate_fractal_motives(bit_len_fractal_colab, 2)}")
            print(f"  Number of motives of size 3: {calculate_fractal_motives(bit_len_fractal_colab, 3)}")

            print("\n--- Demonstrating Law V: AGC-256 Cube Derivation ---")
            print("Simulating AGC-256 Cube derivation from a core nucleotide and context bits:")
            derived_cube_colab = simulate_agc256_cube_derivation(core_nucleotide='A', left_context_bit='0', right_context_bit='1')
            print(derived_cube_colab)


        except Exception as e:
            print(f"Error during V2 Unicode example: {e}")

        # V1 ASCII Encoding Example (for comparison, only works for ASCII characters)
        print("\n---\n")
        ascii_sample_text = "Hello, Colab!"
        print(f"Original ASCII Text (V1 ASCII): {ascii_sample_text}")
        try:
            encoded_v1 = string_to_nucleotide_sequence_v1(ascii_sample_text)
            print(f"Encoded (V1 ASCII): {''.join(encoded_v1)}")

            encoded_v1_with_checksum = add_genetic_checksum(encoded_v1)
            print(f"Encoded with Checksum (V1 ASCII): {''.join(encoded_v1_with_checksum)}")
            print(f"Checksum for V1 is valid: {verify_genetic_checksum(encoded_v1_with_checksum)}")

            decoded_v1 = decode_nucleotide_sequence_to_string_v1(encoded_v1)
            print(f"Decoded (V1 ASCII): {decoded_v1}")
            print(f"V1 Encoding/Decoding successful: {ascii_sample_text == decoded_v1}")
        except Exception as e:
            print(f"Error during V1 ASCII example: {e}")

    else:
        try:
            setup_gui()
        except tk.TclError as e:
            print(f"Error: {e}")
            print("Tkinter GUI cannot be displayed in this environment (e.g., Google Colab). Not a local environment.")
            print("Run this script locally on your computer with a graphical interface.")
