"""The evaluation profile the validator suite runs under: MTB4K_TEST_PROFILE (default r8).

Tests only. report_contract.py activates the profile its process names with --profile, once, at
import, and a test process's own arguments name none, so every test file imports this module
before any lane module that loads report_contract. It imports report_contract with
`--profile PROFILE` as the process's arguments and puts the real arguments back. Under r8 (the
default) nothing changes. Under rN the suite runs against revision_n.py's filled values and the
tree's own checker blobs, as on the revision-9 fallback merge; an unfilled rN refuses here, at
import, as it does in every CLI (tests/test_evaluation_profiles.py checks that with the switch
set). argv(...) adds the same --profile to an in-process CLI call or a subprocess, whose
argparse default is r8.
"""
import os
import sys
from pathlib import Path

LANE = Path(__file__).resolve().parents[1]
if str(LANE) not in sys.path:
    sys.path.insert(0, str(LANE))
PROFILE = os.environ.get("MTB4K_TEST_PROFILE", "r8")
_saved = sys.argv[:]
sys.argv = [_saved[0] if _saved else "test", "--profile", PROFILE]
try:
    import report_contract as vc  # noqa: E402
finally:
    sys.argv = _saved
if vc.PROFILE != PROFILE:
    raise SystemExit(f"profile_boot: report_contract was loaded under {vc.PROFILE}, not MTB4K_TEST_PROFILE={PROFILE}; "
                     "import profile_boot before any lane module")


def exec_module(spec, module):
    """spec.loader.exec_module(module) for a copy of report_contract.py (a mutant), whose import-time
    activation then reads this suite's --profile, not the test process's own arguments."""
    saved = sys.argv[:]
    sys.argv = [saved[0] if saved else "test", "--profile", PROFILE]
    try:
        spec.loader.exec_module(module)
    finally:
        sys.argv = saved


def argv(args):
    """A CLI argument list with this suite's --profile appended."""
    return [*map(str, args), "--profile", PROFILE]
