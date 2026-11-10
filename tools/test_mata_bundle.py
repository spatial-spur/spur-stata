"""Check that the installed Mata bundle is reproducible and self-contained."""
from pathlib import Path
import re
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


class MataBundleTests(unittest.TestCase):
    def test_generation_is_reproducible(self):
        bundle = ROOT / "spur_mata_source.do"
        before = bundle.read_bytes()
        subprocess.run([sys.executable, str(ROOT / "tools/build_mata_bundle.py")], check=True)
        self.assertEqual(before, bundle.read_bytes())

    def test_installer_selects_source_bundle(self):
        manifest = (ROOT / "spur.pkg").read_text()
        entries = [line.split(maxsplit=1)[1] for line in manifest.splitlines() if line.startswith("f ")]
        self.assertNotIn("lspur.mlib", entries)
        for path in entries:
            self.assertTrue((ROOT / path).is_file(), path)
        self.assertIn("spur_mata_source.do", entries)
        self.assertIn("spur_mata.ado", entries)

    def test_public_commands_use_bundled_namespaces(self):
        bundle = (ROOT / "spur_mata_source.do").read_text()
        for path in ROOT.glob("*.ado"):
            if path.name == "spur_mata.ado":
                continue
            for name in re.findall(r"\b(sv2_\w+)\s*\(", path.read_text()):
                self.assertLessEqual(len(name), 32)
                self.assertRegex(bundle, r"\b" + name + r"\s*\(")
        for name in re.findall(r"^struct\s+(sv2_\w+)\s*\{", bundle, re.M):
            self.assertLessEqual(len(name), 32)


if __name__ == "__main__":
    unittest.main()
