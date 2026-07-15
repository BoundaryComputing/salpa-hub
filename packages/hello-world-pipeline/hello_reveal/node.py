"""
hello-reveal — Generates a visual comparison report for the Secret Message pipeline.

Part of the "Secret Message" hello-world workflow template.
Reads both encrypted and decrypted files, then produces a formatted
visual report showing the cipher in action.

All parameters are independently configurable.
Only case_name is received from the upstream node.

No external dependencies — uses Python stdlib only.
Runs IN_PROCESS (no pixi.toml needed).
"""

import os
from datetime import datetime

from bocoflow_core.node import Node, NodeException, NodeResult
from bocoflow_core.parameters import (
    FolderParameter,
    SelectParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log


def build_box_report(case_name, ciphertext, plaintext, shift, timestamp):
    """Build a box-drawing visual report."""
    content_lines = [
        f"Case:       {case_name}",
        f"Cipher:     Caesar (shift={shift})",
        f"Encrypted:  {ciphertext}",
        f"Decrypted:  {plaintext}",
        f"Match:      {'YES' if len(ciphertext) == len(plaintext) else 'length mismatch'}",
    ]
    if timestamp:
        content_lines.append(f"Time:       {timestamp}")

    max_width = max(len(line) for line in content_lines)
    box_width = max(max_width + 4, 42)

    lines = []
    title = " SECRET MESSAGE REVEALED "
    pad = box_width - len(title) - 2
    left_pad = pad // 2
    right_pad = pad - left_pad

    lines.append("+" + "=" * box_width + "+")
    lines.append("|" + " " * left_pad + title + " " * right_pad + "|")
    lines.append("+" + "-" * box_width + "+")

    for content_line in content_lines:
        padding = box_width - len(content_line) - 2
        lines.append("|  " + content_line + " " * padding + "|")

    lines.append("+" + "=" * box_width + "+")

    return "\n".join(lines)


def build_detailed_report(case_name, ciphertext, plaintext, shift, timestamp):
    """Build a detailed report with character mapping and statistics."""
    lines = []

    lines.append("=" * 50)
    lines.append("  SECRET MESSAGE REVEALED")
    lines.append("=" * 50)
    lines.append("")
    lines.append(f"  Case:    {case_name}")
    lines.append(f"  Cipher:  Caesar (shift={shift})")
    if timestamp:
        lines.append(f"  Time:    {timestamp}")
    lines.append("")
    lines.append("-" * 50)
    lines.append(f"  Encrypted:  {ciphertext}")
    lines.append(f"  Decrypted:  {plaintext}")
    lines.append("")

    # Character mapping
    lines.append("-" * 50)
    lines.append("  CHARACTER MAPPING:")
    lines.append("")
    display_len = min(len(ciphertext), 40)
    enc_chars = ciphertext[:display_len]
    dec_chars = plaintext[:display_len]
    lines.append(f"    Encrypted:  {enc_chars}")
    arrows = "".join("|" if c.isalpha() else " " for c in enc_chars)
    lines.append(f"                {arrows}")
    lines.append(f"    Decrypted:  {dec_chars}")
    lines.append("")

    # Statistics
    total_chars = len(plaintext)
    alpha_chars = sum(1 for c in plaintext if c.isalpha())
    shifted_chars = sum(1 for a, b in zip(ciphertext, plaintext) if a != b)

    lines.append("-" * 50)
    lines.append("  STATISTICS:")
    lines.append(f"    Total characters:   {total_chars}")
    lines.append(f"    Alphabetic:         {alpha_chars}")
    lines.append(f"    Shifted:            {shifted_chars}")
    lines.append(f"    Unchanged:          {total_chars - shifted_chars}")
    lines.append("")

    # Round-trip verification
    re_encrypted = []
    for ch in plaintext:
        if ch.isalpha():
            base = ord("A") if ch.isupper() else ord("a")
            re_encrypted.append(chr((ord(ch) - base + shift) % 26 + base))
        else:
            re_encrypted.append(ch)
    match = "".join(re_encrypted) == ciphertext

    lines.append(f"  Verification: {'PASS' if match else 'FAIL'}")
    lines.append("=" * 50)

    return "\n".join(lines)


class HelloReveal(Node):
    """
    Generates a visual comparison report for encrypted vs decrypted messages.

    Input: Predecessor data with case_name (to locate the encrypted/decrypted files)
    Output: A formatted visual report file
    """

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            docstring="Name for this case (used to find the encrypted/decrypted files). "
            "If left empty, uses predecessor data or defaults to 'secret'.",
        ),
        "input_dir": FolderParameter(
            "Input Directory",
            docstring="Directory where the encrypted and decrypted files are located.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            docstring="Directory where the reveal report will be written.",
        ),
        "report_style": SelectParameter(
            "Report Style",
            options=["box", "detailed"],
            default="detailed",
            docstring="Visual style for the report: "
            "'box' for a compact framed summary, "
            "'detailed' for character-by-character analysis with statistics.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting reveal report...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()

            # Read parameters
            case_name = flow_vars["case_name"].get_value()
            input_dir = self.resolve_path(flow_vars["input_dir"].get_value())
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())
            report_style = flow_vars["report_style"].get_value() or "detailed"
            shift = 3  # Default; could be made configurable if needed

            # Use predecessor case_name if ours is empty
            if not case_name and predecessor_data:
                input_data = predecessor_data[0] if predecessor_data else None
                if input_data and isinstance(input_data, dict):
                    case_name = input_data.get("case_name", "secret")
            if not case_name:
                case_name = "secret"

            result.metadata.update({
                "case_name": case_name,
                "execution_time": datetime.now().isoformat(),
            })

            stream_log(f"Reading files for case '{case_name}'...", node_id=self.node_id, progress=20)

            # Read encrypted file (plain ciphertext)
            encrypted_file = os.path.join(input_dir, f"{case_name}_encrypted.txt")
            if not os.path.exists(encrypted_file):
                raise FileNotFoundError(
                    f"Encrypted file not found: {encrypted_file}\n"
                    f"Make sure hello-encrypt ran first with case_name='{case_name}'."
                )

            with open(encrypted_file, "r") as f:
                ciphertext = f.read().strip()

            # Read decrypted file (plain plaintext)
            decrypted_file = os.path.join(input_dir, f"{case_name}_decrypted.txt")
            if not os.path.exists(decrypted_file):
                raise FileNotFoundError(
                    f"Decrypted file not found: {decrypted_file}\n"
                    f"Make sure hello-decrypt ran first with case_name='{case_name}'."
                )

            with open(decrypted_file, "r") as f:
                plaintext = f.read().strip()

            if not ciphertext:
                raise ValueError("Encrypted file is empty.")
            if not plaintext:
                raise ValueError("Decrypted file is empty.")

            stream_log(f"Building {report_style} report...", node_id=self.node_id, progress=50)

            # Generate the report
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

            if report_style == "box":
                report = build_box_report(case_name, ciphertext, plaintext, shift, timestamp)
            else:
                report = build_detailed_report(case_name, ciphertext, plaintext, shift, timestamp)

            # Write output file
            os.makedirs(output_dir, exist_ok=True)
            output_file = os.path.join(output_dir, f"{case_name}_reveal.txt")
            with open(output_file, "w") as f:
                f.write(report + "\n")

            stream_log(f"Wrote: {case_name}_reveal.txt", node_id=self.node_id, progress=90)

            # Build result
            output_path = self.format_output_path(output_file)
            result.data = {
                "case_name": case_name,
                "output_file": output_path,
            }

            result.files["output"] = {
                "report": output_path,
            }

            result.success = True
            result.message = f"Reveal report generated for {case_name}"

            stream_log("Reveal complete", node_id=self.node_id, progress=100)
            return result.to_json()

        except Exception as e:
            stream_log(f"Error: {str(e)}", node_id=self.node_id, progress=100)
            raise NodeException("hello-reveal", str(e))
