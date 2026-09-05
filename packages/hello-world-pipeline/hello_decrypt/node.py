"""
hello-decrypt — Decrypts a Caesar-cipher message from a file.

Part of the "Secret Message" hello-world workflow template.
Reads the encrypted file produced by hello-encrypt, decrypts it,
and writes the plaintext to a new file.

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
    IntegerParameter,
    StringParameter,
)
from bocoflow_core.stream_logger import stream_log


def caesar_decrypt(text, shift):
    """Decrypt text using a Caesar cipher with the given shift."""
    result = []
    for char in text:
        if char.isalpha():
            base = ord("A") if char.isupper() else ord("a")
            result.append(chr((ord(char) - base - shift) % 26 + base))
        else:
            result.append(char)
    return "".join(result)


#: How to run this node on its own -- the values `salpa smoke` feeds it. Strings
#: starting with `demo_data/` resolve relative to this directory; a parameter's
#: type gives its shape and never its value, so every value is written down here.
DEMO_CONFIG = {
    "case_name": 'demo',
    "shift": 3,
    "input_dir": 'demo_data/demo',
}


class HelloDecrypt(Node):
    """
    Decrypts a Caesar-cipher message from an encrypted file.

    Input: Predecessor data with case_name (to locate the encrypted file)
    Output: A text file containing the decrypted plaintext
    """

    OPTIONS = {
        "case_name": StringParameter(
            "Case Name",
            docstring="Name for this case (used to find the encrypted file). "
            "If left empty, uses predecessor data or defaults to 'secret'.",
        ),
        "shift": IntegerParameter(
            "Cipher Shift",
            default=3,
            docstring="Number of positions to reverse-shift each letter (1-25). "
            "Must match the shift used for encryption.",
        ),
        "input_dir": FolderParameter(
            "Input Directory",
            default="rel:./",
            docstring="Directory where the encrypted file is located.",
        ),
        "output_dir": FolderParameter(
            "Output Directory",
            default="rel:./",
            docstring="Directory where the decrypted file will be written.",
        ),
    }

    def execute(self, predecessor_data, flow_vars):
        stream_log("Starting decryption...", node_id=self.node_id, progress=0)

        try:
            result = NodeResult()

            # Read parameters
            case_name = flow_vars["case_name"].get_value()
            shift = flow_vars["shift"].get_value()
            if shift is None:
                shift = 3
            shift = max(1, min(25, int(shift)))
            input_dir = self.resolve_path(flow_vars["input_dir"].get_value())
            output_dir = self.resolve_path(flow_vars["output_dir"].get_value())

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

            stream_log(f"Reading {case_name}_encrypted.txt...", node_id=self.node_id, progress=20)

            # Read the encrypted file (plain ciphertext, one line)
            encrypted_file = os.path.join(input_dir, f"{case_name}_encrypted.txt")
            if not os.path.exists(encrypted_file):
                raise FileNotFoundError(
                    f"Encrypted file not found: {encrypted_file}\n"
                    f"Make sure hello-encrypt ran first with case_name='{case_name}'."
                )

            with open(encrypted_file, "r") as f:
                ciphertext = f.read().strip()

            if not ciphertext:
                raise ValueError("Encrypted file is empty.")

            stream_log(f"Decrypting with shift={shift}...", node_id=self.node_id, progress=50)

            plaintext = caesar_decrypt(ciphertext, shift)

            stream_log(f"'{ciphertext}' -> '{plaintext}'", node_id=self.node_id, progress=70)

            # Write plain decrypted text file
            os.makedirs(output_dir, exist_ok=True)

            output_file = os.path.join(output_dir, f"{case_name}_decrypted.txt")
            with open(output_file, "w") as f:
                f.write(plaintext + "\n")

            stream_log(f"Wrote: {case_name}_decrypted.txt", node_id=self.node_id, progress=90)

            # Build result
            output_path = self.format_output_path(output_file)
            result.data = {
                "case_name": case_name,
                "output_file": output_path,
            }

            result.files["output"] = {
                "decrypted": output_path,
            }

            result.success = True
            result.message = f"Message decrypted for {case_name}"

            stream_log("Decryption complete", node_id=self.node_id, progress=100)
            return result.to_json()

        except Exception as e:
            stream_log(f"Error: {str(e)}", node_id=self.node_id, progress=100)
            raise NodeException("hello-decrypt", str(e))
