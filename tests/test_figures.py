from pathlib import Path
from microdim_repro.cli import reference_figures


def test_reference_figures(tmp_path):
    reference_figures(tmp_path)
    expected=['Figure_1.png','Figure_2.png','Figure_3.png','Figure_4.png','Figure_S1.png','Figure_S2.png','Figure_S3.png','Figure_S4.png']
    for x in expected:
        p=tmp_path/x
        assert p.exists() and p.stat().st_size>1000
