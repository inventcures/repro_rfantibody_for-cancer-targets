"""Patch RFAntibody's get_next_frames to handle degenerate rotation matrices.

RFdiffusion can predict zero rotation matrices for certain target structures,
causing scipy.spatial.transform.Rotation.from_matrix to crash with
'ValueError: Non-positive determinant'. This patch replaces degenerate matrices
(det < 1e-6) with identity matrices before passing to scipy.

Patches both R_0 and R_t calls in get_next_frames().
"""

import pathlib
import re
import sys

UTILS_PATH = pathlib.Path(
    "/opt/RFAntibody/src/rfantibody/rfdiffusion/inference/utils.py"
)

TARGETS = [
    ("R_0", "R_0 = scipy_R.from_matrix(R_0.squeeze().numpy()).as_matrix()"),
    ("R_t", "R_t = scipy_R.from_matrix(R_t.squeeze().numpy()).as_matrix()"),
]


def _make_replacement(var: str, indent: str) -> str:
    tmp = f"{var}_np"
    return "\n".join([
        f"{indent}{tmp} = {var}.squeeze().numpy()",
        f"{indent}dets = np.linalg.det({tmp})",
        f"{indent}bad = np.abs(dets) < 1e-6",
        f"{indent}if bad.any():",
        f"{indent}    {tmp}[bad] = np.eye(3)",
        f"{indent}{var} = scipy_R.from_matrix({tmp}).as_matrix()",
    ])


def patch():
    if not UTILS_PATH.exists():
        print(f"ERROR: {UTILS_PATH} not found", file=sys.stderr)
        sys.exit(1)

    text = UTILS_PATH.read_text()
    patched = 0

    for var, old_line in TARGETS:
        if old_line not in text:
            print(f"  {var}: already patched or not found — skipping")
            continue

        match = re.search(r"^( +)" + re.escape(old_line), text, re.MULTILINE)
        if not match:
            print(f"  {var}: could not detect indentation — skipping")
            continue

        indent = match.group(1)
        text = text.replace(match.group(0), _make_replacement(var, indent), 1)
        patched += 1
        print(f"  {var}: patched ({len(indent)}-space indent)")

    if patched:
        UTILS_PATH.write_text(text)
        print(f"Patched {patched} scipy_R.from_matrix calls in get_next_frames")
    else:
        print("Nothing to patch")


if __name__ == "__main__":
    patch()
