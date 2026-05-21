import subprocess
import sys
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "verify_code.py"


class VerifyCodeCliTests(unittest.TestCase):
    def run_cmd(self, *args, cwd=None):
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), *args],
            cwd=str(cwd or REPO_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        return proc.returncode, proc.stdout + proc.stderr

    def test_default_target_passes(self):
        code, out = self.run_cmd()
        self.assertEqual(code, 0, msg=out)
        self.assertIn("=== 结果: ✅ 全部通过 ===", out)

    def test_non_r_target_is_rejected(self):
        code, out = self.run_cmd("--file", "verify_code.py")
        self.assertEqual(code, 3, msg=out)
        self.assertIn("不支持的文件类型", out)

    def test_missing_target_returns_1(self):
        code, out = self.run_cmd("--file", "does_not_exist.R")
        self.assertEqual(code, 1, msg=out)
        self.assertIn("目标文件不存在", out)

    def test_subdir_invocation_still_works(self):
        code, out = self.run_cmd(cwd=REPO_ROOT / "modules")
        self.assertEqual(code, 0, msg=out)
        self.assertIn("目标文件:", out)


if __name__ == "__main__":
    unittest.main()
