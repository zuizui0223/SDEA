"""Generate SDEA v18 retrospective evidence ladder from earlier source-locked receipts."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sdea.evidence_ladder import main

if __name__ == "__main__":
    main()
