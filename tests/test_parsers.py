import os
import tempfile
from siren_rnai.siren_plotIV import parse_fasta

def test_parse_fasta():
    """Test that the FASTA parser correctly extracts sequences and lengths."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        tmp.write(">Mock_Gene_01\nATGCGTACGT\n")
        tmp_path = tmp.name
    
    try:
        sequence, length = parse_fasta(tmp_path)
        assert sequence == "ATGCGTACGT"
        assert length == 10
    finally:
        os.remove(tmp_path)
