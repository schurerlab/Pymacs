from __future__ import annotations

import ast
import itertools
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path("/Users/jxs794/Documents/Pymacs")
TARGET_FILES = [
    ROOT / "3A_AutomateGromacs.py",
    ROOT / "3A_AutomateGromacs_MPI.py",
    ROOT / "3A_AutomateGromacs_chunks.py",
    ROOT / "3A_AutomateGromacs_chunks_MPI.py",
]

FUNCTION_NAMES = [
    "normalize_chain_type",
    "is_polymer_chain_type",
    "chain_display_label",
    "collect_polymer_chain_entries",
    "_interface_chain_pairs_for_manifest",
]


def load_functions(script_path: Path):
    tree = ast.parse(script_path.read_text(encoding="utf-8"), filename=str(script_path))
    selected = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in FUNCTION_NAMES:
            selected.append(node)
    module = ast.Module(body=selected, type_ignores=[])
    namespace = {"combinations": itertools.combinations}
    exec(compile(module, str(script_path), "exec"), namespace)
    return namespace


class InterfaceDisplayLabelTests(unittest.TestCase):
    def test_legacy_entries_without_display_label_build_manifest_pairs(self):
        chain_map = {"1kgy": (0, 1000), "xEFN3_2dock": (1001, 2000)}
        chain_entry_meta = {
            "1kgy": {
                "label": "1kgy",
                "start": 0,
                "end": 1000,
                "current_chain": None,
                "chain_type": "protein",
            },
            "xEFN3_2dock": {
                "label": "xEFN3_2dock",
                "start": 1001,
                "end": 2000,
                "current_chain": None,
                "chain_type": "protein",
            },
        }
        for path in TARGET_FILES:
            with self.subTest(script=path.name):
                ns = load_functions(path)
                entries = ns["collect_polymer_chain_entries"](chain_map, chain_entry_meta, {})
                self.assertEqual(entries[0]["display_label"], "1kgy")
                self.assertEqual(entries[1]["display_label"], "xEFN3_2dock")
                ns["polymer_chain_entries"] = entries
                ns["args"] = SimpleNamespace(interface_chain_pairs=None)
                self.assertEqual(ns["_interface_chain_pairs_for_manifest"](), ["1kgy:xEFN3_2dock"])

    def test_current_chain_takes_precedence(self):
        entry = {"label": "Protein1", "current_chain": "A", "display_label": "OldLabel"}
        for path in TARGET_FILES:
            with self.subTest(script=path.name):
                ns = load_functions(path)
                self.assertEqual(ns["chain_display_label"](entry, 1), "A")

    def test_source_chain_fallback(self):
        entry = {"label": "Protein1", "current_chain": None, "source_chain": "SRC"}
        for path in TARGET_FILES:
            with self.subTest(script=path.name):
                ns = load_functions(path)
                self.assertEqual(ns["chain_display_label"](entry, 1), "SRC")

    def test_existing_display_label_beats_plain_label(self):
        entry = {"label": "foo", "display_label": "nice_label"}
        for path in TARGET_FILES:
            with self.subTest(script=path.name):
                ns = load_functions(path)
                self.assertEqual(ns["chain_display_label"](entry, 1), "nice_label")

    def test_sparse_entry_still_has_label(self):
        entry = {"label": "foo"}
        for path in TARGET_FILES:
            with self.subTest(script=path.name):
                ns = load_functions(path)
                self.assertEqual(ns["chain_display_label"](entry, 3), "foo")

    def test_final_fallback_never_crashes(self):
        entry = {}
        for path in TARGET_FILES:
            with self.subTest(script=path.name):
                ns = load_functions(path)
                self.assertEqual(ns["chain_display_label"](entry, 7), "Chain_7")


if __name__ == "__main__":
    unittest.main()
